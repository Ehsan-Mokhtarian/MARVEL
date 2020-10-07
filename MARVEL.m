function [G_out, mb_tests, marvel_tests] = MARVEL(varargin)

% MARVEL recovers a dag in the Markov equivalence class of the 
% ground-truth dag. 
% Use the function in the following way with name-value pairs.
% [G, mb_tests, marvel_tests] = MARVEL('data', D, 'alpha',0.01,...
% 'initial_mb','TC', 'check_Mb_inclusion',false, ...
% 'top_down_subsets',true, 'use_L0', true, ... 
% 'cond_indep','stat', 'G',zeros(p))

% Input arguments: 
% 'data': n*p data matrix - n samples of each of the p variables. necessary
% if the CI tests are performed on data. (not necessary for oracle CI
% tests)
% alpha: the parameter for fischer Z-transform - default value:0.01
% initial_mb: the algorithm to recover the Markov blanket information. 
% possible values: 'TC' (default), 'GS', 'TCbw'.
% check_Mb_inclusion: true-false value indicating whether to do the 
% preprocess of checking if \Mb_X \subseteq \Mb_Y for every Y\in\Mb_X.
% possible values: true (default), false
% top_down_subsets:  true-false value indicating whether to start testing
% the subsets of a set from small to large subsets or vice versa.
% possible values: true, false (default)
% use_L0: true-false value indicating whether or not to save the result of 
% the performed CI tests to avoid duplicate tests. possible values: 
% true (default), false
% cond_indep: Type of the performed CI tests. use 'oracle' for oracle CI 
% tests. in this case you have to feed the ground-truth graph G. use 
% 'stat' (default) for statistical CI tests, which is Fischer's Z 
% transform. (CI tests are implemented for Gaussian variables for now.)
% G: the adjacency matrix of the ground truth. Necessary in the oracle 
% setting to perform the oracle CI tests.
% Output values:
% G_out: a dag recovered by MARVEL, which is in the MEC of the ground truth
% mb_tests: number of CI tests required to compute the Markov blankets
% marvel_tests: number of CI tests that MARVEL performs given Markov 
% blanket information

%--------------------------------------------------------------------------
%************************ Default Config **********************************

p = 10;
n = 100;
MARVEL_Options = struct('data', zeros(n,p),'alpha',0.01,'initial_mb','TC',...
    'check_Mb_inclusion',false,'top_down_subsets', true,...
    'use_L0', true, 'cond_indep', 'stat', 'G', zeros(p));
optionNames = fieldnames(MARVEL_Options);
% Finalize Config according to the input arguments:
for option_pair = reshape(varargin,2,[])
   if any(strcmp(option_pair{1},optionNames))
      MARVEL_Options.(option_pair{1}) = option_pair{2};
   else
      error('Unrecognized parameter %s',option_pair{1})
   end
end
D = MARVEL_Options.('data');
alpha = MARVEL_Options.('alpha');
initial_mb = MARVEL_Options.('initial_mb');
check_Mb_inclusion = MARVEL_Options.('check_Mb_inclusion');
top_down_subsets = MARVEL_Options.('top_down_subsets');
use_L0 = MARVEL_Options.('use_L0');
cond_indep = MARVEL_Options.('cond_indep');
G = MARVEL_Options.('G');
if strcmp(cond_indep, 'oracle')
    p = size(G,2);
    D = zeros(n,p);
else
    p = size(D,2);
    G = zeros(p);
end

%*********************** Start of the algorithm ***************************

V = ones(1,p);
E = zeros(p);
L_0 = zeros(p);
L_2 = zeros(p,p,p);
separating_sets = cell(1,0);
extra_links = false(p);
marvel_tests = cell(1,0);
[Mb, mb_tests] = Markov_Blanket(G,D, initial_mb, cond_indep, 2/p/(p-1));
mb_tests = length(mb_tests);
% Markov blanket information discovered


for iter=p:-1:1
    mbs = sum(Mb);
    [~,ind] = sort(mbs);
    XX = 0;
    ind = ind(V(ind)>0);
    thresh = min(p,length(ind));
    counter = thresh;
    version = 1;
    if counter == 0 
        version = 2;
    end
    f1 = false(thresh,1);
    for X=ind
        if version ==2 
            [out_f1, Mb, L_0, extra_links, separating_sets,marvel_tests] = Filter1(G,D, X, Mb, check_Mb_inclusion, top_down_subsets, use_L0, L_0, extra_links, separating_sets,alpha,marvel_tests,cond_indep);
            if out_f1
                XX = X;
                [Mb, separating_sets, extra_links, L_2,marvel_tests] = Filter2_ver2(G,D,X,Mb,top_down_subsets,separating_sets,extra_links,use_L0,L_2,alpha,marvel_tests,cond_indep);
                break
            end 
        else
            [f1(thresh+1-counter), Mb, L_0, extra_links, separating_sets,marvel_tests] = Filter1(G,D, X, Mb, check_Mb_inclusion, top_down_subsets, use_L0, L_0, extra_links, separating_sets,alpha,marvel_tests,cond_indep);
            if f1(thresh+1-counter) % X passes Filter1
                [f2, Mb, L_0, L_2, extra_links, separating_sets,marvel_tests] = Filter2(G,D, X, Mb, top_down_subsets, L_0, L_2, use_L0, extra_links, separating_sets,alpha,marvel_tests,cond_indep);
                if f2 % X passes Filter2
                    XX = X;
                    break
                end
            end
        end
        if XX == 0 % no removable found
            counter = counter - 1;
            if counter == 0
                version = 2;
                elligibles = find(f1);
                if ~isempty(elligibles)
                    XX = ind(elligibles(1));
                    [Mb, separating_sets, extra_links, L_2,marvel_tests] = Filter2_ver2(G,D,XX,Mb,top_down_subsets,separating_sets,extra_links,use_L0,L_2,alpha,marvel_tests,cond_indep);
                    break
                end
            end
        end
    end
    if XX == 0 
        XX = ind(1);
        [Mb, extra_links, separating_sets, L_0,marvel_tests] = Filter1_ver2(G,D,XX,Mb,top_down_subsets,extra_links,separating_sets,use_L0,L_0,alpha,marvel_tests,cond_indep);
        [Mb, separating_sets, extra_links, L_2,marvel_tests] = Filter2_ver2(G,D,XX,Mb,top_down_subsets,separating_sets,extra_links,use_L0,L_2,alpha,marvel_tests,cond_indep);
    end
    if version == 1
        [Mb,E,V,L_0,marvel_tests] = Update_Graph(G,D, XX, Mb, E, V, L_0,alpha,marvel_tests,cond_indep);
    else
        [Mb,E,V,L_0,marvel_tests] = Update_Graph_ver2(G,D,XX,Mb, E, V, L_0,alpha,marvel_tests,cond_indep);
    end
end

G_out = (E - extra_links)>0;

% Correct the v-structures if there are any extra_links: (modified version of MARVEL)

for i=1:length(separating_sets)
    set = separating_sets{i};
    X = set(1);
    Y = set(2);
    comm_neighbors = myintersect(myunion(find(G_out(X,:)),find(G_out(:,X))), myunion(find(G_out(Y,:)),find(G_out(:,Y))));
    sep_set = mysetdiff(set,[X Y]);
    for W=comm_neighbors
        if ~ismember(W, sep_set)
            G_out(W,X) = false;
            G_out(W,Y) = false;
            G_out(X,W) = true;
            G_out(Y,W) = true;
        end
    end
end
marvel_tests = length(marvel_tests);
end
%************************* End of Algorithm *******************************
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%**************************** Functions ***********************************


function [out, Mb, L_0, extra_links, separating_sets, tests] = Filter1(G,D, X, Mb, check_Mb_inclusion, top_down_subsets, use_L0, L_0, extra_links, separating_sets,alpha,tests,cond_indep)
    out = true;
    mb = find(Mb(X,:));
    if check_Mb_inclusion
        for Y=mb
            if ~all(ismember([mb X], [find(Mb(Y,:)) Y]))  %Mb(X) subseteq Mb(Y)
                out = false;
                return
            end
        end
    end
    
    if use_L0
        if any(L_0(X,:)>0)
            out = false;
            return
        end
    end
    
    mbs = size(mb,2);
    for Yind=1:mbs
        Y=mb(Yind);
        if use_L0 && L_0(X,Y) == -1 && L_0(Y,X) == -1
            continue
        end
        AA =[mb(1:(Yind-1)), mb((Yind+1):mbs)];
        CI = false;
        if top_down_subsets
            subset_sizes = (mbs-2):-1:0;
        else
            subset_sizes = 0:(mbs-2);
        end
        for sizes = subset_sizes
            mb_subsets = subsets1(AA,sizes);
            for t=1:length(mb_subsets)
                S = mb_subsets{t};
                if strcmp(cond_indep, 'oracle')
                    CI = CI_Test(X,Y,S,D,alpha,cond_indep,G);
                else
                    CI = CI_Test(X,Y,S,D,alpha,cond_indep);
                end
                tests{end+1} = [min(X,Y) max(X,Y) S];
                if CI
                    out = false;
                    L_0(X,Y) = 1;
                    L_0(Y,X) = 1;
                    separating_sets{end+1} = [X, Y, S]; 
                    extra_links(Y,X) = true;
                    extra_links(X,Y) = true;
                    return
                end
            end
        end
        
        if ~CI
            L_0(X,Y) = -1;
            L_0(Y,X) = -1;
        end
    end
end

function [out, Mb, L_0, L_2, extra_links, separating_sets,tests] = Filter2(G,D, X, Mb, top_down_subsets, L_0, L_2, use_L0, extra_links, separating_sets,alpha,tests,cond_indep)
    out = true;
    mb = find(Mb(X,:));
    mbs = size(mb,2);
    if use_L0 && any(any(L_2(X,mb,mb)))
        out=false;
        return
    end
    
    for Yind=1:(mbs-1)    
        Y=mb(Yind);
        for Wind=(Yind+1):mbs    
            W=mb(Wind);
            CI = false;
            if use_L0 && L_2(X,Y,W) == -1 && L_2(X,W,Y) == -1
                continue
            end
            AA =[mb(1:(Yind-1)) ,mb((Yind+1):(Wind-1)), mb((Wind+1):mbs)];    
            if top_down_subsets
                subset_sizes = (length(AA)):-1:0;
            else
                subset_sizes = 0:(length(AA));
            end
            for sizes = subset_sizes
                mb_subsets = subsets1(AA,sizes);
                for t=1:length(mb_subsets)
                    S = mb_subsets{t};
                    if strcmp(cond_indep, 'oracle')
                        CI = CI_Test(Y,W,[X,S],D,alpha,cond_indep,G);
                    else
                        CI = CI_Test(Y,W,[X,S],D,alpha,cond_indep);
                    end
                    tests{end+1} = [min(Y,W) max(Y,W) [X,S]];
                    if CI       
                        L_2(X,Y,W) = 1;
                        L_2(X,W,Y) = 1;
                        if Mb(Y,W) && Mb(W,Y)
                            L_0(Y,W) = 1;
                            L_0(W,Y) = 1;
                            separating_sets{end+1} = [Y, W, [X,S]]; 
                            extra_links(Y,W) = true;
                            extra_links(W,Y) = true;
                        end
                        out = false;
                        return;
                    end
                end
            end
            
            
            if ~CI
                L_2(X,Y,W) = -1;
                L_2(X,W,Y) = -1;
            end
        end
    end
end

function [Mb,E,V,L_0,tests] = Update_Graph(G,D, X, Mb, E, V, L_0,alpha,tests,cond_indep)
    mb = find(Mb(X,:));
    mbs = size(mb,2);

    % Add New Edges
    for Y=mb
        E(Y,X)=1;
    end

    % Update V
    V(X) = 0;

    % Update Markov Blankets        
    for Y=mb
        Mb(X,Y) = 0;
        Mb(Y,X) = 0;
        L_0(X,Y) = 0;
        L_0(Y,X) = 0;
    end
    for Yind=1:(mbs-1)
        Y = mb(Yind);
        for Wind=(Yind+1):mbs
            W = mb(Wind);
            if nnz(Mb(Y,:))>nnz(Mb(W,:))
                S = mysetdiff(find(Mb(W,:)),Y);
            else
                S = mysetdiff(find(Mb(Y,:)),W);
            end
            if strcmp(cond_indep, 'oracle')
                CI = CI_Test(Y,W,S,D,alpha,cond_indep,G);
            else
                CI = CI_Test(Y,W,S,D,alpha,cond_indep);
            end
            tests{end+1} = [min(Y,W) max(Y,W) S];
            if CI
                Mb(W,Y) = 0;
                Mb(Y,W) = 0;
                L_0(W,Y) = 0;
                L_0(Y,W) = 0;
            end
        end
    end
end

function [Mb, separating_sets, extra_links, L_2,tests] = Filter2_ver2(G,D,X,Mb,top_down_subsets,separating_sets,extra_links,use_L0,L_2,alpha,tests,cond_indep)
    mb = find(Mb(X,:));
    mbs = size(mb,2);
    for Yind=1:(mbs-1)    
        Y=mb(Yind);
        for Wind=(Yind+1):mbs    
            W=mb(Wind);     
            if use_L0 && L_2(X,Y,W) == -1 && L_2(X,W,Y) == -1
                continue
            end
            if extra_links(Y,W)
                continue
            end
            AA =[mb(1:(Yind-1)) ,mb((Yind+1):(Wind-1)), mb((Wind+1):mbs)];    
            
            if top_down_subsets
                subset_sizes = (mbs-2):-1:0;
            else
                subset_sizes = 0:(mbs-2);
            end
            for sizes = subset_sizes
                mb_subsets = subsets1(AA,sizes);
                for t=1:length(mb_subsets)
                    S = mb_subsets{t};
                    if strcmp(cond_indep, 'oracle')
                        CI = CI_Test(Y,W,[X,S],D,alpha,cond_indep,G);
                    else
                        CI = CI_Test(Y,W,[X,S],D,alpha,cond_indep);
                    end
                    tests{end+1} = [min(Y,W) max(Y,W) [X,S]];
                    if CI
                        separating_sets{end+1} = [Y, W, [X,S]]; 
                        extra_links(Y,W) = true;
                        extra_links(W,Y) = true;
                        break;
                    end
                end
            end
            
            if CI
                continue;
            end
        end
    end
end

function [Mb,E,V,L_0,tests] = Update_Graph_ver2(G,D,X,Mb, E, V, L_0,alpha,tests,cond_indep)
    mb = find(Mb(X,:));
    mbs = size(mb,2);

    % Add New Edges
    for Y=mb
        E(Y,X)=1;
    end

    % Update V
    V(X) = 0;

    % Update Markov Blankets        
    for Y=mb
        Mb(X,Y) = 0;
        Mb(Y,X) = 0;
        L_0(X,Y) = 0;
        L_0(Y,X) = 0;
    end
    for Yind=1:(mbs-1)
        Y = mb(Yind);
        for Wind=(Yind+1):mbs
            W = mb(Wind);
            if L_0(W,Y) == 1 && L_0(Y,W) == 1
                L_0(W,Y) = 0;
                L_0(Y,W) = 0;
            end
            if nnz(Mb(Y,:))>nnz(Mb(W,:))
                S = mysetdiff(find(Mb(W,:)),Y);
            else
                S = mysetdiff(find(Mb(Y,:)),W);
            end
            if strcmp(cond_indep, 'oracle')
                CI = CI_Test(Y,W,S,D,alpha,cond_indep,G);
            else
                CI = CI_Test(Y,W,S,D,alpha,cond_indep);
            end
            tests{end+1} = [min(Y,W) max(Y,W) S];
            if CI
                Mb(W,Y) = 0;
                Mb(Y,W) = 0;
                L_0(W,Y) = 0;
                L_0(Y,W) = 0;
            end
        end
    end
end

function [Mb, extra_links, separating_sets, L_0,tests] = Filter1_ver2(G,D,X,Mb,top_down_subsets,extra_links,separating_sets,use_L0,L_0,alpha,tests,cond_indep)
    mb = find(Mb(X,:));
    mbs = size(mb,2);
    
    for Yind=1:mbs
        Y=mb(Yind);
        if use_L0 && L_0(X,Y) == -1 && L_0(Y,X) == -1
            continue
        end
        if extra_links(X,Y)
            continue
        end
        AA =[mb(1:(Yind-1)), mb((Yind+1):mbs)];
        
        if top_down_subsets
            subset_sizes = (mbs-2):-1:0;
        else
            subset_sizes = 0:(mbs-2);
        end
        for sizes = subset_sizes
            mb_subsets = subsets1(AA,sizes);
            for t=1:length(mb_subsets)
                S = mb_subsets{t};
                if strcmp(cond_indep, 'oracle')
                    CI = CI_Test(X,Y,S,D,alpha,cond_indep,G);
                else
                    CI = CI_Test(X,Y,S,D,alpha,cond_indep);
                end
                tests{end+1} = [min(X,Y) max(X,Y) S];
                if CI
                    separating_sets{end+1} = [X, Y, S]; 
                    extra_links(Y,X) = true;
                    extra_links(X,Y) = true;
                    break
                end
            end
        end
    end
end

function [Mb, tests] = Markov_Blanket(G,D, initial_mb, cond_indep, alpha)
    tests = cell(1,0); % CI tests done
    if strcmp(cond_indep, 'oracle')
        p = size(G,2);
    else
        p = size(D,2);
    end
    Mb = zeros(p);
    if strcmp(initial_mb, 'TC') % Mb Through Total Conditioning
        [Mb,tests] = total_conditioning(G,D, alpha,tests,cond_indep);
    elseif strcmp(initial_mb, 'GS') % GS algo
        [Mb,tests] = grow_shrink(G,D, alpha,tests,cond_indep);
    elseif strcmp(initial_mb, 'TCbw') % Mb Through Backward Total Conditioning
        Mb = TCbw(D, 0.05, alpha);
    end
end

function [Mb,tests] = total_conditioning(G,D, alpha,tests,cond_indep)
% Total Conditioning, see Pellet and Elisseeff
    if strcmp(cond_indep, 'oracle')
        p = size(G,2);
    else
        p = size(D,2);
    end
    Mb = zeros(p);
    for X=1:(p-1)
        for Y=(X+1):p
            S = [1:(X-1),(X+1):(Y-1),(Y+1):p];
            if strcmp(cond_indep, 'oracle')
                CI = CI_Test(X,Y,S,D,alpha,cond_indep,G);
            else
                CI = CI_Test(X,Y,S,D,alpha,cond_indep);
            end
            tests{end+1} = [min(X,Y) max(X,Y) S];
            if ~CI
                Mb(X,Y)=1;
                Mb(Y,X)=1;
            end
        end
    end
end


function [Mb,tests] = grow_shrink(G,D, alpha,tests,cond_indep)
% GS algorithm. See Margaritis and Thrun
    if strcmp(cond_indep, 'oracle')
        p = size(G,2);
    else
        p = size(D,2);
    end
    Mb = zeros(p);
    for X=1:p
        has_changed = true;
        while has_changed
            has_changed = false;
            other_nodes = mysetdiff([1:(X-1),(X+1):p], find(Mb(X,:)));
            for Y = other_nodes  % Grow Phase
                S = find(Mb(X,:));
                if strcmp(cond_indep, 'oracle')
                    CI = CI_Test(X,Y,S,D,alpha,cond_indep,G);
                else
                    CI = CI_Test(X,Y,S,D,alpha,cond_indep);
                end
                tests{end+1} = [min(X,Y) max(X,Y) S];
                if ~CI
                    Mb(X,Y)=1;
                    has_changed = true;
                end
            end
        end
        has_changed = true;
        while has_changed
            has_changed = false;
            mb = find(Mb(X,:));
            for Y=mb  % Shrink Phase
                S = mysetdiff(mb, Y);
                if strcmp(cond_indep, 'oracle')
                    CI = CI_Test(X,Y,S,D,alpha,cond_indep,G);
                else
                    CI = CI_Test(X,Y,S,D,alpha,cond_indep);
                end
                tests{end+1} = [min(X,Y) max(X,Y) S];
                if  CI
                    Mb(X,Y) = 0;
                    has_changed = true;
                end
            end
        end
    end
    % One last symmetry check:
    for X=1:p
        for Y=(X+1):p
            if ~Mb(X,Y) || ~Mb(Y,X)
                Mb(X,Y) = 0;
                Mb(Y,X) = 0;
            end
        end
    end
end


function Mb = TCbw(threshold)
% backward TC. See Pellet and Elisseeff
    global p D
    Mb = zeros(p);
    for X = 1:p
        P = [1:(X-1) (X+1):p];
        while ~isempty(P) && length(P)>nnz(Mb(X,:))
            weights = regress(D(:,X), D(:,P));
            [~,ind] = sort(weights,'descend');
            least_significant_pred = P(ind(ceil(length(P)/2):length(P)));
            sig_pred = P(weights>threshold);
            Mb(X,sig_pred) = 1;
            P = mysetdiff(P, least_significant_pred);
        end
    end
    % One last symmetry check:
    for X=1:p
        for Y=(X+1):p
            if ~Mb(X,Y) || ~Mb(Y,X)
                Mb(X,Y) = 0;
                Mb(Y,X) = 0;
            end
        end
    end
end


% From here on, we have used a few basic functions written by other people,
% namely Kevin Murphy et al. See
% https://github.com/bayesnet/bnt/find/master.


function C = mysetdiff(A,B)
% MYSETDIFF Set difference of two sets of positive integers (much faster than built-in setdiff)
% C = mysetdiff(A,B)
% C = A \ B = { things in A that are not in B }
%
% Original by Kevin Murphy, modified by Leon Peshkin

if isempty(A)
    C = [];
    return;
elseif isempty(B)
    C = A;
    return; 
else % both non-empty
    bits = zeros(1, max(max(A), max(B)));
    bits(A) = 1;
    bits(B) = 0;
    C = A(logical(bits(A)));
end
end


function C = myunion(A,B)
% MYUNION Union of two sets of positive integers (much faster than built-in union)
% C = myunion(A,B)

if isempty(A)
  ma = 0;
else
  ma = max(A);
end

if isempty(B)
  mb = 0;
else
  mb = max(B);
end

if ma==0 && mb==0
  C = [];
elseif ma==0 && mb>0
  C = B;
elseif ma>0 && mb==0
  C = A;
else
  %bits = sparse(1, max(ma,mb));
  bits = zeros(1, max(ma,mb));
  bits(A) = 1;
  bits(B) = 1;
  C = find(bits);
end
end


function C = myintersect(A,B)
% MYINTERSECT Intersection of two sets of positive integers (much faster than built-in intersect)
% C = myintersect(A,B)

A = A(:)'; B = B(:)';

if isempty(A)
  ma = 0;
else
  ma = max(A);
end

if isempty(B)
  mb = 0;
else
  mb = max(B);
end

if ma==0 || mb==0
  C = [];
else
  %bits = sparse(1, max(ma,mb));
  bits = zeros(1, max(ma,mb));
  bits(A) = 1;
  C = B(logical(bits(B)));  
end
end


function sub_s=subsets1(s,k)
% SUBSETS1 creates sub-sets of a specific from a given set
% SS = subsets1(S, k)
% 
% S is the given set
% k is the required sub-sets size
% 
% Example:
% 
% >> ss=subsets1([1:4],3);
% >> ss{:}
% ans =
%      1     2     3
% ans =
%      1     2     4
% ans =
%      1     3     4
% ans =
%      2     3     4
% 
% Written by Raanan Yehezkel, 2004

if k<0 % special case
    error('subset size must be positive');
elseif k==0 % special case
    sub_s={[]};
else
    l=length(s);
    ss={};
    if l>=k
        if k==1 % Exit condition
            for I=1:l
                ss{I}=s(I);
            end
        else
            for I=1:l
                ss1=subsets1(s([(I+1):l]),k-1);
                for J=1:length(ss1)
                    ss{end+1}=[s(I),ss1{J}];
                end
            end
        end
    end
    sub_s=ss;
end
end

%***************************** The End ************************************
%--------------------------------------------------------------------------