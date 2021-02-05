function [A,D,perm_matrix] = Generate_Data(varargin)


    DATA_options = struct('graph_type', 'mildew', 'sample_size', 500,...
        'num_of_variables',25, 'max_indegree',5,'density',0.1, ...
        'num_of_edges',15);
    optionNames = fieldnames(DATA_options);
    % Finalize Config according to the input arguments:
    for option_pair = reshape(varargin,2,[])
       if any(strcmp(option_pair{1},optionNames))
          DATA_options.(option_pair{1}) = option_pair{2};
       else
          error('Unrecognized parameter %s',option_pair{1})
       end
    end
    graph_type = DATA_options.('graph_type');
    n = DATA_options.('sample_size');
    p = DATA_options.('num_of_variables'); 
    d = DATA_options.('max_indegree');
    r = DATA_options.('density');
    edges = DATA_options.('num_of_edges');

    [A, perm_matrix] = Generate_Graph(graph_type,p,d,r,edges);
    p = size(A,1);
    N = randn(n,p)*diag(1 + 2*rand(1,p));
    AA = (0.5+0.5*rand(p)).* ((-1).^(rand(p)>0.5));
    AA = A.*AA;
    D = N/(eye(p)-AA);
end



