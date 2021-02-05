function [A, perm_matrix] = Generate_Graph(graph_type,p,d,r,edges)
    %Generate_Graph generates the adjacency matrix of a DAG. if graph type is
    %random, the output graoh will have p vertices and max indegree d.
    % 'random', 'diabetes', 'alarm', 'insurance', 'hailfinder', 'carpo'
    if strcmp(graph_type, 'random_num_of_edges_fixed')
        AA = triu(rand(p),1);
        [~,sortingIndex] = sort(AA(:), 'desc');
        [row, col] = ind2sub(size(AA),sortingIndex(1:edges));
        A = zeros(p);
        for i=1:length(row)
            A(row(i),col(i)) = 1;
        end
    elseif strcmp(graph_type, 'random_sparsity_fixed')
        A = rand(p)<r;
        A = triu(A,+1);
    elseif strcmp(graph_type, 'random_maxindegree_fixed')
        A = zeros(p);
        for i=2:p
            x = randperm(p)';           
            A(x(1:d),i) = 1;         
        end
        A = triu(A,+1);    
    else
        f = matfile(graph_type+".mat");
        A = f.A;
    end
    % randomly permute A:
    p = size(A, 1);
    perm_matrix = eye(p);
    perm_matrix = perm_matrix(randperm(p),:);
    A = perm_matrix*A*(perm_matrix');
end
