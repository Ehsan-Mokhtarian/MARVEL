function sep = dsep(G, X, Y, S)
% DSEP Is X indep Y given S wrt DAG G?
% sep = dsep(X, Y, S, G)
%
% Instead of using the Bayes-Ball criterion, we see if S separates X and Y
% in the moralized ancestral graph.
conn = reachability_graph(G);
M = myunion(myunion(X, Y), S);
[AA,~] = find(conn(:, M));
AA = unique(AA);
AA = myunion(AA, M);
GM = moralize(G(AA,AA));
%sep = graph_separated(GM, X, Y, S);
sep = graph_separated(GM, find_equiv_posns(X,AA), find_equiv_posns(Y,AA), find_equiv_posns(S,AA));
end


function C = reachability_graph(G)
% REACHABILITY_GRAPH C(i,j) = 1 iff there is a path from i to j in DAG G
% C = reachability_graph(G)

if 1
    % expm(G) = I + G + G^2 / 2! + G^3 / 3! + ...
    M = expm(double(full(G))) - eye(length(G));
    C = (M>0.00001);
else
    % This computes C = G + G^2 + ... + G^{n-1}
    n = length(G);
    A = G;
    C = zeros(n);
    for i=1:n-1
        C = C + A;
        A = A * G;
    end
    C = (C > 0);
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


function [M, moral_edges] = moralize(G)
% MORALIZE Ensure that for every child, all its parents are married, and drop directionality of edges.
% [M, moral_edges] = moralize(G)

M = G;
n = length(M);
for i=1:n
    fam = family(G,i);
    M(fam,fam)=1;
end
M = setdiag(M,0);
moral_edges = sparse(triu(max(0,M-G),1));
end

function f = family(A,i,t)
% FAMILY Return the indices of parents and self in sorted order
% f = family(dag,i,t)
%
% t is an optional argument: if present, dag is assumed to be a 2-slice DBN

if nargin < 3
    f = [parents(A,i) i];
else
    if t == 1
        f = [parents(A,i) i];
    else
        ss = length(A)/2;
        j = i+ss;
        f = [parents(A,j) j] + (t-2)*ss;
    end
end
end

function ps = parents(adj_mat, i)
% PARENTS Return the list of parents of node i
% ps = parents(adj_mat, i)

ps = find(adj_mat(:,i))';
end


function M = setdiag(M, v)
% SETDIAG Set the diagonal of a matrix to a specified scalar/vector.
% M = set_diag(M, v)

n = length(M);
if length(v)==1
    v = repmat(v, 1, n);
end

% e.g., for 3x3 matrix,  elements are numbered
% 1 4 7
% 2 5 8
% 3 6 9
% so diagnoal = [1 5 9]


J = 1:n+1:n^2;
M(J) = v;

%M = triu(M,1) + tril(M,-1) + diag(v);
end


function p = find_equiv_posns(vsmall, vlarge)
% FIND_EQUIV_POSNS p[i] = the place where vsmall[i] occurs in vlarge.
% p = find_equiv_posns(vsmall, vlarge)
% THE VECTORS ARE ASSUMED TO BE SORTED.
%
% e.g., vsmall=[2,8], vlarge=[2,7,8,4], p=[1,3]
%
% In R/S, this function is called 'match'

%if ~mysubset(vsmall, vlarge)
%  error('small domain must occur in large domain');
%end

if isempty(vsmall) || isempty(vlarge)
    p = [];
    return;
end

bitvec = sparse(1, max(vlarge));
%bitvec = zeros(1, max(vlarge));
bitvec(vsmall) = 1;
p = find(bitvec(vlarge));

%p = find(ismember(vlarge, vsmall)); % slower
end

function sep = graph_separated(G, X, Y, S)

G2 = G;
G2(S,:) = 0;
G2(:,S) = 0;
conn = reachability_graph(G2);
conn2 = conn(X,Y);
sep = all(conn2(:)==0);
end

