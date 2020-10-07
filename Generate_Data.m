function [A,D,perm_matrix] = Generate_Data(graph_type,n)
    [A, perm_matrix] = Generate_Graph(graph_type);
    p = size(A,1);
    N = randn(n,p)*diag(1 + 2*rand(1,p));
    AA = ((-1).^(rand(p)>0.5)).*(0.5+1.5*rand(p));
    AA = A.*AA;
    D = N/(eye(p)-AA);
end


function [A, perm_matrix] = Generate_Graph(graph_type,p,d,r,edges)
    %Generate_Graph generates the adjacency matrix of a DAG. if graph type is
    %random, the output graoh will have p vertices and max indegree d.
    % 'random', 'diabetes', 'alarm', 'insurance', 'hailfinder', 'carpo'
    if nargin < 2
        d=1;
        r=0.1;
        edges = 10;
        p = 10;
    end
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
%         A = zeros(p);
%         for i=2:p
%             x = rand(i-1,1)<(d/(p-1));
%             while sum(x)>d
%                 x = rand(i-1,1)<(d/(p-1));
%             end
%             A(1:(i-1),i) = x;
%         end
    elseif strcmp(graph_type, 'alarm')
        f = matfile('alarm.mat');
        A = f.A;
    elseif strcmp(graph_type, 'andes')
        f = matfile('andes.mat');
        A = f.A;
    elseif strcmp(graph_type, 'asia')
        f = matfile('asia.mat');
        A = f.A;
    elseif strcmp(graph_type, 'barley')
        f = matfile('barley.mat');
        A = f.A;
    elseif strcmp(graph_type, 'cancer')
        f = matfile('cancer.mat');
        A = f.A;
    elseif strcmp(graph_type, 'carpo')
        f = matfile('carpo.mat');
        A = f.A;
    elseif strcmp(graph_type, 'diabetes')
        f = matfile('diabetes.mat');
        A = f.A;
    elseif strcmp(graph_type, 'hailfinder')
        f = matfile('hailfinder.mat');
        A = f.A;
    elseif strcmp(graph_type, 'hepar2')
        f = matfile('hepar2.mat');
        A = f.A;
    elseif strcmp(graph_type, 'insurance')
        f = matfile('insurance.mat');
        A = f.A;
    elseif strcmp(graph_type, 'Link')
        f = matfile('Link.mat');
        A = f.A;
    elseif strcmp(graph_type, 'mehra')
        f = matfile('mehra.mat');
        A = f.A;
	elseif strcmp(graph_type, 'mildew')
        f = matfile('mildew.mat');
        A = f.A;
    elseif strcmp(graph_type, 'munin1')
        f = matfile('munin1.mat');
        A = f.A;
    elseif strcmp(graph_type, 'munin2')
        f = matfile('munin2.mat');
        A = f.A;
    elseif strcmp(graph_type, 'munin3')
        f = matfile('munin3.mat');
        A = f.A;
    elseif strcmp(graph_type, 'munin4')
        f = matfile('munin4.mat');
        A = f.A;
    elseif strcmp(graph_type, 'pathfinder')
        f = matfile('pathfinder.mat');
        A = f.A;
    elseif strcmp(graph_type, 'pigs')
        f = matfile('pigs.mat');
        A = f.A;
    elseif strcmp(graph_type, 'sachs')
        f = matfile('sachs.mat');
        A = f.A;
    elseif strcmp(graph_type, 'water')
        f = matfile('water.mat');
        A = f.A;
    elseif strcmp(graph_type, 'win95pts')
        f = matfile('win95pts.mat');
        A = f.A;
    else
        A = zeros(p);
    end
    % randomly permute A:
    p = size(A, 1);
    perm_matrix = eye(p);
    perm_matrix = perm_matrix(randperm(p),:);
    A = perm_matrix*A*(perm_matrix');
end

