function [extra_edges, missing_edges, precision, recall, skeleton_F1_score] = learning_errors(A,E_algo, perm_matrix)
    
    p = size(A, 1);
    if nargin<1
        E_algo = zeros(p);
    end
    
    A = standard_pdag(perm_matrix'*A*perm_matrix);
    E_algo = standard_pdag(perm_matrix'*E_algo*perm_matrix);

    [extra_edges, missing_edges, precision, recall] = skeleton_errors(A, E_algo);
    skeleton_F1_score = 2*precision*recall/(precision+recall);
end


function [extra_edges, missing_edges, precision, recall] = skeleton_errors(G, H_algo)
    % G is the ground truth pdag, while H_algo is the pdag an algo
    % recovers.
    G = abs(G);
    G = (G+G')>0;   % build undirected graphs
    H_algo = abs(H_algo);
    H_algo = (H_algo+H_algo')>0;
    skeleton_errors = H_algo - G;
    extra_edges = nnz(skeleton_errors>0)/2;
    missing_edges = nnz(skeleton_errors<0)/2;
    precision = nnz(G.*H_algo)/nnz(H_algo);
    recall = nnz(G.*H_algo)/nnz(G);
end



function pdag = standard_pdag(G)
%standard_pdag converts a pdag to its standard format:
%   pdag(i,j) = -1 if there is an i->j edge.
%   pdag(i,j) = G(j,i) = 1 if there is an undirected edge i <-> j
%   input G MUST be in the right topoligal order
    p = size(G,1);
    G = abs(G);
    G = (G+G')>0;   % undirected graph (skeleton)
    pdag = ones(p).*G;
    for X = 3:p
        parents = find(G(X,1:X));
        if length(parents)<2
            continue
        end
        couples = nchoosek(parents,2);
        for i=1:size(couples,1)
            P1 = couples(i,1);
            P2 = couples(i,2);
            if ~G(P1,P2)
                pdag(P1,X) = -1;
                pdag(P2,X) = -1;
                pdag(X,P1) = 0;
                pdag(X,P2) = 0;
            end
        end
    end
end

