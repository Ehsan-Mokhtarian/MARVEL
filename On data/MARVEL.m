function [G_MARVEL,tests,SC] = MARVEL(D,Mb,alpha)

% MARVEL algorithm
% input:
% D: data matrix. Each colomn corresponds to a variable, and each row
% corresponds to a sample
% Mb: The matrix of Markov boundaries
% alpha: significance level for CI tests
%
% Output:
% G_MARVEL: The learned graph using MARVEL algorithm
% tests: The number of performed CI tests
% SC: A vector of length p+1 such that the i-th entry indicates the
% number of performed CI tests with a conditioning set of size i-1

%*********************** Initialization ***************************
p = size(D,2);
M1 = -ones(p).*Mb; % Learned skeleton, 1: edge, 0: no edge, -1: unknown
M2 = false(p,p,p); % Z = M2(:,X,Y) where X d-sep Y|Z
M3 = zeros(p,p,p); % M3(X,T,Z)=1 if T and Z does not have any separating set in Mb_X \cup X.

V = ones(1,p);
G_MARVEL = zeros(p);
tests = 0;
SC = zeros(p+1,1);

%*********************** Start of the algorithm ***************************
for iter=p:-1:1
    Mb_size = sum(Mb);
    [~,ind] = sort(Mb_size);
    XX = 0; % the removable node that we want to find
    ind = ind(V(ind)>0);
    
    % Finding XX
    for X=ind
        Mb_X = find(Mb(:,X));
        [isR,M1,M2,M3,N,VS,tests,SC] = IsRemovable...
            (X,Mb_X,G_MARVEL,D,alpha,M1,M2,M3,tests,SC);
        if isR
            XX = X;
            break
        end
    end
    
    % Updating\Removing XX
    if XX == 0
        %         disp(remaining_nodes)
        G_MARVEL(ind, ind) = M1(ind, ind);
        break
    else
        [Mb,G_MARVEL,V,M1,M2,tests,SC] = UpdateGraph...
            (XX,M1,M2,N,VS,D,Mb,G_MARVEL,V,alpha,tests,SC);
    end
end
end

%************************* End of Algorithm *******************************
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%**************************** Functions ***********************************

function [isR,M1,M2,M3,N,VS,tests,SC] = IsRemovable...
    (X,Mb_X,E,D,alpha,M1,M2,M3,tests,SC)

VS = [];  % V-Structures
isR = true;
p = size(D,2);

% Condition 1:
[is_neighbors,sep_sets,M1,M2,tests,SC] = neighbors...
    (X,Mb_X,D,alpha,M1,M2,tests,SC);
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
                CI = CI_Test(W,Y,[S X],D,0.7);
                tests = tests+1;
                ls = length(S)+2;
                SC(ls) =  SC(ls)+1;
                if CI
                    isR = false;
                    sep = false(p,1);
                    sep([S X])=true;
                    M2(:,W,Y) = sep;
                    M2(:,Y,W) = sep;
                    M1(W, Y) = 0;
                    M1(Y, W) = 0;
                    return
                end
            end
        end
        M3(X, W, Y) = 1;
        M3(X, Y, W) = 1;
    end
end

% Condition 2
[VS, M1,M2,M3,tests,SC] = v_structures...
    (X,Mb_X,E,is_neighbors,sep_sets,D,alpha,M1,M2,M3,tests,SC);
M4 = zeros(p,p,p); % M4(T,Z,Y)=1 if T and Z do not have any separating set including Y in Mb_X \cup X
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
        for sep_size = 0:length(search_set)
            mb_subsets = subsets1(search_set,sep_size);
            for t=1:length(mb_subsets)
                S = mb_subsets{t};
                CI = CI_Test(W,T,[S X Y],D,0.7);
                tests = tests+1;
                ls = length(S)+3;
                SC(ls) = SC(ls)+1;
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


function [VS,M1,M2,M3,tests,SC] = v_structures...
    (X,Mb_X,E,is_neighbors,sep_sets,D,alpha,M1,M2,M3,tests,SC)
%v_structures finds the set of v_structures of a given vertex X.

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
        [is_vstructs(i,j),M1,M2,M3,tests,SC] = IsVstructure...
            (X,neighbor,cp,S,Mb_X,D,alpha,M1,M2,M3,tests,SC);
    end
end
VS = zeros(sum(is_vstructs,'all'), 2);
[indices_i, indices_j] = find(is_vstructs);
for k = 1:length(indices_i)
    VS(k,1) = Mb_X(neighbors_X(indices_i(k)));
    VS(k,2) = Mb_X(coparents_X(indices_j(k)));
end
end

function [is_vstructure,M1,M2,M3,tests,SC] = IsVstructure...
    (X,Z,T,S,Mb_X,D,alpha,M1,M2,M3,tests,SC)
% is X - Z - T a vstructure? given that Z is a neighbor of X and T is a
% coparent of X. S is the separating set for X,T.

p = size(D,2);
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
    CI = CI_Test(Z,T,S,D,alpha);
    tests = tests+1;
    ls = length(S)+1;
    SC(ls) = SC(ls)+1;
    if CI
        is_vstructure = false;
        M1(T,Z) = 0;
        M1(Z,T) = 0;
        sep = false(p,1);
        sep(S) = true;
        M2(:,Z,T) = sep;
        M2(:,T,Z) = sep;
        return;
    end
    
    search_set = mysetdiff(Mb_X, [T, Z]);
    for sep_size = 0:length(search_set)
        mb_subsets = subsets1(search_set,sep_size);
        for t=1:length(mb_subsets)
            sep_set = mb_subsets{t};
            CI = CI_Test(Z,T,[sep_set X],D,alpha);
            ls = length(sep_set)+2;
            SC(ls) = SC(ls)+1;
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


function [is_neighbors,sep_sets,M1,M2,tests,SC] = neighbors...
    (X,Mb_X,D,alpha,M1,M2,tests,SC)
%neighbors finds the set of neighbors of a given vertex X.

p = size(D,2);
mbs = length(Mb_X);
is_neighbors = false(1, mbs);
sep_sets = false(p,mbs);
for i = 1:length(Mb_X)
    T = Mb_X(i);
    [is_neighbors(i), sep_sets(:,i), M1,M2,tests,SC] = IsNeighbor...
        (X,T,Mb_X,D,alpha,M1,M2,tests,SC);
end
end

function [is_neighbor,sep_set,M1,M2,tests,SC] = IsNeighbor...
    (X,T,Mb_X,D,alpha,M1,M2,tests,SC)

%is_neighbor determines whether X,T are neighbors. If not, it also returns
%sep_set \subset mb_X that separates X,T.

p = size(D,2);
if M1(X, T) == 0  % We already know that they are not neighbors!
    Z = M2(:,X,T);
    if is_subset(find(Z),Mb_X)
        is_neighbor = false;
        sep_set = Z;
        return
    end
elseif M1(X, T) == 1  % We already know that they are neighbors!
    is_neighbor = true;
    sep_set = false(p,1);
    return
end
search_set = mysetdiff(Mb_X, T);
for sep_size = 0:(length(search_set)-1)
    mb_subsets = subsets1(search_set,sep_size);
    for t=1:length(mb_subsets)
        S = mb_subsets{t};
        CI = CI_Test(X,T,S,D,alpha);
        tests = tests + 1;
        ls = length(S)+1;
        SC(ls) = SC(ls)+1;
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


function [Mb,E,V,M1,M2,tests,SC] = UpdateGraph...
    (XX,M1,M2,N,VS,D,Mb,E,V,alpha,tests,SC)

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
    return
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
        CI = CI_Test(Y,W,S,D,alpha);
        ls = length(S)+1;
        SC(ls) = SC(ls)+1;
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