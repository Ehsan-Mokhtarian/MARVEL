function [E, Mb_tests, MARVEL_tests] = MARVEL(varargin)

% MARVEL recovers a dag in the Markov equivalence class of the 
% ground-truth dag.
% Use the function in the following way with name-value pairs.
% [G, mb_tests, marvel_tests] = MARVEL('data', D, 'alpha',0.01,...
% 'Mb_algorithm','TC','cond_indep','stat', 'G',zeros(p))

%************************* Input arguments*********************************
% 'data': n*p data matrix - n samples, p variables. necessary if the CI 
% tests are performed on data. (not necessary for oracle CI tests)

% 'alpha': the parameter for fischer Z-transform - default value:0.01

% 'Mb_algorithm': the algorithm to recover the Markov blanket information. 
% possible values: 'TC' (default), 'GS'.

% 'cond_indep': Type of the performed CI tests. use 'oracle' for oracle CI 
% tests. in this case you have to feed the ground-truth graph G. use 'stat'
% (default) for statistical CI tests, which is Fischer's Z transform.
% (CI tests are implemented for Gaussian variables for now.)

% 'G': the adjacency matrix of the ground truth. Necessary in the oracle 
% setting to perform the oracle CI tests.

% 'plot_on': Plot the reamining graph at each step. Only for 'oracle'.
% true or false (Default = false)

%************************* Output arguments********************************
% E: Adjacency matrix of a dag recovered by MARVEL, which is in the MEC of 
% the ground truth

% Mb_tests: number of CI tests required to compute the Markov blankets

% MARVEL_tests: number of CI tests that MARVEL performs given Markov 
% blanket information

%************************* Default Config *********************************
% Default options
p = 10;
n = 100;
MARVEL_Options = struct('data', zeros(n,p),'alpha',0.01,...
    'Mb_algorithm','TC','cond_indep','stat','plot_on',false,'G',zeros(p));
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
Mb_algorithm = MARVEL_Options.('Mb_algorithm');
cond_indep = MARVEL_Options.('cond_indep');
G = MARVEL_Options.('G');
plot_on = MARVEL_Options.('plot_on');
if strcmp(cond_indep, 'oracle')
    p = size(G,2);
    D = zeros(n,p);
else
    p = size(D,2);
    G = zeros(p);
end

%*********************** Start of the algorithm ***************************

M1 = -ones(p); % Learned skeleton, 1: edge, 0: no edge, -1: unknown
M2 = false(p,p,p); % Z = M2(:,X,Y) where X d-sep Y|Z
M3 = zeros(p,p,p); % M3(X,T,Z)=1 if T and Z does not have any 
                   % separating set in Mb_X \cup X.

alpha_Mb = 2/(p^2);
[Mb, Mb_tests] = Markov_Blanket(D, Mb_algorithm, alpha_Mb, cond_indep,G);
M1 = M1.*Mb;  

V = ones(1,p);
E = zeros(p);
MARVEL_tests = 0;

for iter=p:-1:1
    Mb_size = sum(Mb);
    [~,ind] = sort(Mb_size);
    XX = 0; % the removable node that we want to find
    ind = ind(V(ind)>0);

    % Finding XX
    for X=ind
        Mb_X = find(Mb(:,X));
        [isR,M1,M2,M3,N,VS,MARVEL_tests] = IsRemovable(X,...
            Mb_X,E,D,alpha,cond_indep,M1,M2,M3,MARVEL_tests,G);
        if isR
            XX = X;                
            break
        end 
    end

    remaining_nodes = find(V);
    if plot_on
        if iter==p
            figure
            P1 = plot(digraph(G(V==1,V==1),string(remaining_nodes)));
            xdata = P1.XData;
            ydata = P1.YData;
            XL = xlim;
            YL = ylim;
        else
            P1 = plot(digraph(G(V==1,V==1),string(remaining_nodes)), ...
                'XData',xdata(V==1), 'YData',ydata(V==1));
            xlim(XL);
            ylim(YL);
        end       
        highlight(P1,find(remaining_nodes==XX,1),'NodeColor','r')
        highlight(P1,find(remaining_nodes==XX,1))
        pause
    end


    % Updating\Removing XX
    if XX == 0
        disp(remaining_nodes)
        E(ind, ind) = M1(ind, ind);
        break;
    else
        alpha_update = alpha_Mb;%2/(sum(V)^2);
        [Mb,E,V,M1,M2,MARVEL_tests] = UpdateGraph...
            (XX,M1,M2,N,VS, G, D, Mb, E, V,alpha_update,...
            MARVEL_tests,cond_indep);
    end
end  
end

%************************* End of Algorithm *******************************
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%**************************** Functions ***********************************

function [isR,M1,M2,M3,N,VS,tests] = IsRemovable...
    (X,Mb_X,E,D,alpha,cond_indep,M1,M2,M3,tests,varargin)
    
    VS = [];  % V-Structures
    isR=true;
    alpha_removable = 0.8;
    if strcmp(cond_indep, 'oracle')
        G = varargin{1};
        p = size(G,2);
    else
        G = 0;
        p = size(D,2);
    end

    % Condition 1:    
    [is_neighbors, sep_sets, M1,M2,tests] = neighbors(X, Mb_X, D, ...
        alpha, cond_indep, M1, M2,tests, G);
    N = Mb_X(is_neighbors); % Neighbors
    NS = length(N); % Neighbor size
    for i=1:(NS-1)
        W = N(i);
        for j=(i+1):NS
            Y = N(j);
            if M1(W,Y)==1 || M3(X, W, Y) == 1
                continue
            end
            search_set = mysetdiff(Mb_X, [W Y]);
            for sep_size = 0:length(search_set)
                mb_subsets = subsets1(search_set,sep_size);
                for t=1:length(mb_subsets)
                    S = mb_subsets{t};
                    CI = CI_Test(W,Y,[S X],D,alpha_removable,cond_indep,G);
                    tests = tests+1;
                    if CI
                        isR = false;
                        sep = false(p,1);
                        sep([S X])=true;
                        M2(:,W,Y) = sep;
                        M2(:,Y,W) = sep;
                        M1(W, Y) = 0;
                        M1(Y, W) = 0;
                        return;
                    end
                end
            end
            M3(X, W, Y) = 1;
            M3(X, Y, W) = 1;
        end
    end
    
    % Condition 2    
    [VS, M1,M2,M3,tests] = v_structures(X, Mb_X, E, is_neighbors,...
        sep_sets, D, alpha, cond_indep, M1,M2,M3,tests, G);
    M4 = zeros(p,p,p); % M4(T,Z,Y)=1 if T and Z do not have any 
                       % separating set including Y in Mb_X \cup X 
    for i=1:size(VS,1)
        Y = VS(i,1);
        T = VS(i,2);
        for W = mysetdiff(N,Y)'
            if M1(W,T)==1 || M3(X, W, T) == 1
                continue
            end
            if M1(W,T)==0 &&  is_subset([X,Y] ,find(M2(:,W,T)))
                isR = false;
                return
            end
            search_set = mysetdiff(Mb_X, [T Y W find(M4(T,W,:))']);
            disp(find(M4(T,W,:))');
            for sep_size = 0:length(search_set)
                mb_subsets = subsets1(search_set,sep_size);
                for t=1:length(mb_subsets)
                    S = mb_subsets{t};
                    CI = CI_Test...
                        (W,T,[S X Y],D,alpha_removable,cond_indep,G);
                    tests = tests+1;
                    if CI
                        isR = false;
                        sep = false(p,1);
                        sep([S X Y])=true;
                        M2(:,W,T) = sep;
                        M2(:,T,W) = sep;
                        M1(W, T) = 0;
                        M1(T, W) = 0;
                        return;
                    end
                end
            end
            M4(T,W,Y) = 1;
        end
    end
end


function [Mb, tests] = Markov_Blanket...
    (D, Mb_algorithm, alpha, cond_indep, varargin)
    
    tests = 0; % CI tests done       
    if strcmp(cond_indep, 'oracle')
        G = varargin{1};
        p = size(G,2);
    else
        p = size(D,2);
        G = zeros(p);
    end

    Mb = zeros(p);
    if strcmp(Mb_algorithm, 'TC') % Mb Through Total Conditioning
        [Mb,tests] = total_conditioning(G,D, alpha,tests,cond_indep);
    elseif strcmp(Mb_algorithm, 'GS') % GS algo
        [Mb,tests] = grow_shrink(G,D, alpha,tests,cond_indep);
    end
end

function [Mb,tests] = total_conditioning...
    (G,D, alpha,tests,cond_indep)
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
            tests = tests+1;
            if ~CI
                Mb(X,Y)=1;
                Mb(Y,X)=1;
            end
        end
    end
end

function [Mb,tests] = grow_shrink...
    (G,D, alpha,tests,cond_indep)
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
            % Grow Phase
            for Y = other_nodes  
                S = find(Mb(X,:));
                if strcmp(cond_indep, 'oracle')
                    CI = CI_Test(X,Y,S,D,alpha,cond_indep,G);
                else
                    CI = CI_Test(X,Y,S,D,alpha,cond_indep);
                end
                tests = tests+1;
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
            % Shrink Phase
            for Y=mb  
                S = mysetdiff(mb, Y);
                if strcmp(cond_indep, 'oracle')
                    CI = CI_Test(X,Y,S,D,alpha,cond_indep,G);
                else
                    CI = CI_Test(X,Y,S,D,alpha,cond_indep);
                end
                tests = tests+1;
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


function [VS,M1,M2,M3,tests] = v_structures...
    (X,Mb_X,E,is_neighbors,sep_sets,D,alpha,...
    cond_indep,M1,M2,M3,tests,varargin)
%v_structures finds the set of v_structures of a given vertex X.

    if strcmp(cond_indep, 'oracle')
        G = varargin{1};        
    else
        G = 0;
    end

    neighbors_X = find(is_neighbors);
    coparents_X = find(~is_neighbors);
    
    if length(neighbors_X)==1
        Y = Mb_X(neighbors_X);
        VS = zeros(length(coparents_X),2);
        for j = 1:length(coparents_X)   
            cp = Mb_X(coparents_X(j));
            VS(j,1) = Y;
            VS(j,2) = cp;
            M1(Y,cp) = 1;
            M1(cp,Y) = 1;
        end
        return
    end
    
    is_vstructs = false(length(neighbors_X), length(coparents_X));
    for i = 1:length(neighbors_X)
        neighbor = Mb_X(neighbors_X(i));
        if E(neighbor, X) == 1
            continue
        end
        for j = 1:length(coparents_X)
            cp = Mb_X(coparents_X(j));
            if E(neighbor, cp) == 1
                continue
            end
            S = find(sep_sets(:,coparents_X(j)))';
            [is_vstructs(i,j), M1,M2,M3,tests] = IsVstructure(X,neighbor, ...
                cp,S, Mb_X,D,alpha,cond_indep,M1,M2,M3,tests,G);
        end
    end
    VS = zeros(sum(is_vstructs,'all'), 2);
    [indices_i, indices_j] = find(is_vstructs);
    for k = 1:length(indices_i)
        VS(k,1) = Mb_X(neighbors_X(indices_i(k)));
        VS(k,2) = Mb_X(coparents_X(indices_j(k)));
    end
end

function [is_vstructure, M1,M2,M3,tests] = IsVstructure...
    (X, Z, T, S,Mb_X, D, alpha, cond_indep, M1,M2,M3,tests, varargin)
% is X - Z - T a vstructure? given that Z is a neighbor of X and T is a
% coparent of X. S is the separating set for X,T.
    
    if strcmp(cond_indep, 'oracle')
        p = size(varargin{1},2);
    else
        p = size(D,2);
    end
    if is_subset(Z, S)
        is_vstructure = false;
        return
    end
    if M3(X, T, Z) == 0 && M1(T, Z) < 1
        if M1(T, Z) == 0
            is_vstructure = false;
            return
        end
        
        % Check S before other sets
        if strcmp(cond_indep, 'oracle')
            CI = CI_Test(Z,T,S,D,alpha,cond_indep,varargin{1});
        else
            CI = CI_Test(Z,T,S,D,alpha,cond_indep);
        end
        tests = tests+1;
        if CI
            is_vstructure = false;
            M1(T,Z) = 0;
            M1(Z,T) = 0;         
            sep = false(p,1);
            sep(S)=true;
            M2(:,Z,T) = sep;
            M2(:,T,Z) = sep;
            return;
        end
               
        search_set = mysetdiff(Mb_X, [T, Z]);     
        for sep_size = 0:length(search_set)
            mb_subsets = subsets1(search_set,sep_size);
            for t=1:length(mb_subsets)
                sep_set = mb_subsets{t};
                if strcmp(cond_indep, 'oracle')
                    CI = CI_Test...
                        (Z,T,[sep_set X],D,alpha,cond_indep,varargin{1});
                else
                    CI = CI_Test(Z,T,[sep_set X],D,alpha,cond_indep);
                end
                tests = tests+1;
                if CI
                    is_vstructure = false;
                    M1(T,Z) = 0;
                    M1(Z,T) = 0;         
                    sep = false(p,1);
                    sep(sep_set)=true;
                    M2(:,Z,T) = sep;
                    M2(:,T,Z) = sep;
                    return;
                end
            end
        end
    end
    M3(X, T, Z) = 1;
    M3(X, Z, T) = 1;
    M1(T,Z) = 1;
    M1(Z,T) = 1; 
    is_vstructure = true;
end


function [is_neighbors,sep_sets,M1,M2,tests] = neighbors...
    (X,Mb_X,D,alpha,cond_indep,M1,M2,tests,varargin)
%neighbors finds the set of neighbors of a given vertex X.

    if strcmp(cond_indep, 'oracle')
        G = varargin{1};
        p = size(G,2);
    else
        G = 0;
        p = size(D,2);
    end
    
    mbs = length(Mb_X);
    is_neighbors = false(1, mbs);
    sep_sets = false(p,mbs);
    for i = 1:length(Mb_X)
        T = Mb_X(i);
        [is_neighbors(i), sep_sets(:,i), M1,M2,tests] = IsNeighbor(X, T,...
      Mb_X, D, alpha, cond_indep, M1, M2,tests, G);
    end
end

function [is_neighbor,sep_set,M1,M2,tests] = IsNeighbor...
    (X,T,Mb_X,D,alpha,cond_indep,M1,M2,tests,varargin)

%is_neighbor determines whether X,T are neighbors. If not, it also returns
%sep_set \subset mb_X that separates X,T.
    
    if strcmp(cond_indep, 'oracle')
        G = varargin{1};
        p = size(G,2);
    else
        p = size(D,2);
    end
    
    if M1(X, T) == 0  % We already know that they are not neighbors!
        Z = M2(:,X,T);
        if is_subset(find(Z),Mb_X)
            is_neighbor = false;
            sep_set = Z;
            return;
        end
    elseif M1(X, T) == 1  % We already know that they are neighbors!
        is_neighbor = true;
        sep_set = false(p,1);
        return;
    end
    search_set = mysetdiff(Mb_X, T);
    for sep_size = 0:(length(search_set)-1)
        mb_subsets = subsets1(search_set,sep_size);
        for t=1:length(mb_subsets)
            S = mb_subsets{t};
            if strcmp(cond_indep, 'oracle')
                CI = CI_Test(X,T,S,D,alpha,cond_indep,G);
            else
                CI = CI_Test(X,T,S,D,alpha,cond_indep);
            end
            tests = tests+1;
            if CI
                is_neighbor = false;
                sep_set = false(p,1);
                sep_set(S)=true;
                M2(:,X,T) = sep_set;
                M2(:,T,X) = sep_set;
                M1(X, T) = 0;
                M1(T, X) = 0;
                return;
            end
        end
    end

    M1(X, T) = 1;
    M1(T, X) = 1;
    is_neighbor = true;
    sep_set = false(p,1);
end


function [Mb,E,V,M1,M2,tests] = UpdateGraph...
    (XX,M1,M2,N,VS,G,D,Mb,E,V,alpha,tests,cond_indep)

    % Add New Edges
    children = unique(VS(:,1));
    for child = children
        E(XX, child) = 1;
    end
    for pair = VS'
        E(pair(2),pair(1))=1;
    end
    possible_parents = mysetdiff(N, find(E(XX,:)));
    for parent = possible_parents
        E(parent, XX) = 1;
    end
    
    % Update V
    V(XX) = 0;

    % Update Markov Blankets 
    Mb_X = find(Mb(XX,:));
    for Y=Mb_X
        Mb(XX,Y) = 0;
        Mb(Y,XX) = 0;
    end
    if any(E(XX,V==1))
        return;
    end
    
    ns = length(N);  % size of the possible parents set
    for Yind=1:(ns-1)
        Y = N(Yind);
        for Wind=(Yind+1):ns
            W = N(Wind);
            if M1(W,Y) == 1  % They are neighbors
                continue;
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
            tests = tests+1;
            if CI
                Mb(W,Y) = 0;
                Mb(Y,W) = 0;
                p = size(E,1);
                sep_set = false(p,1);
                sep_set(S)=true;
                M2(:,W,Y) = sep_set;
                M2(:,Y,W) = sep_set;
                M1(Y, W) = 0;
                M1(W, Y) = 0;
            end
        end
    end
end


% From here on, we have used a few basic functions written by other people,
% namely Kevin Murphy et al. See
% https://github.com/bayesnet/bnt/find/master.

function C = mysetdiff(A,B)
    % MYSETDIFF Set difference of two sets of positive integers 
    % (much faster than built-in setdiff)
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

function is_subset = is_subset(A, B)
%is_subset determines if A is a subset of B
    is_subset = all(ismember(A, B));
end

%***************************** The End ************************************
%--------------------------------------------------------------------------