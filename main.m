close all
clear
clc
addpath('./structures');


%%

num_of_samples = 2000;
selected_structure = 6; 

% select ground-truth graph
data_structure  =  ["cancer","asia","sachs","insurance","water"...
    ,"mildew","alarm","barley","hailfinder","carpo","hepar2"...
    ,"win95pts","diabetes","pathfinder","andes","pigs","Link"...
    ,"munin2","munin4","munin1","munin3"];
fprintf('selected structure: %s\n\n', data_structure(selected_structure));

% generate data

% data options:
% 1) 'graph_type': 'random_num_of_edges_fixed', 'random_sparsity_fixed',
% 'random_maxindegree_fixed', or one of the data structures above
% default: 'mildew'

% 2) 'sample_size': positive integer, number of samples of each variable
% default: 500

% 3) 'num_of_variables': order of the graph, used only for random graphs
% default: 25
        
% 4) 'max_indegree': positive integer, maximum indegree of the graph, used 
% only in the 'random_maxindegree_fixed' model
% default: 5

% 5) 'density': real number in the interval [0,1], fraction of the edges 
% present in the graph, used only in the 'random_sparsity_fixed' model
% default: 0.1

% 6) 'num_of_edges': positive integer, number of the edges in the graph,
% used only in the 'random_num_of_edges_fixed' model
% default: 15
    
    
data_options = {'graph_type', data_structure(selected_structure), ...
    'sample_size', num_of_samples};
[G,D,perm_matrix] = Generate_Data(data_options{:});

%% How to call MARVEL

%  MARVEL Options:

% 'cond_indep': CI tests, 'oracle' or 'stat' (Default= 'stat'). 
% 'Mb_algorithm': the algorithm to discover the Markov blanket info. 
% algorithms: 'TC', 'GS' (Default ='TC')

% For 'stat':
% 'alpha': the parameter of Fischer's Z transform (Default =0.01)

% For 'oracle':
% 'G': % Ground truth 
% 'plot_on': Plot the reamining graph at each step. Only for 'oracle'.
% true or false (Default = false)



% two examples:
% 
% 1) oracle setting. Ground truth, G, must be fed to do the CI tests.

MARVEL_Opt1 = {'cond_indep', 'oracle', 'G', G,'plot_on',false};
fprintf('orcale setting...\n');
[E1, Mb_tests, MARVEL_tests] = MARVEL(MARVEL_Opt1{:});
fprintf('MARVEL is done. Number of CI tests: %d \n',Mb_tests+MARVEL_tests);
[extra_edges, missing_edges, precision, recall, skeleton_F1_score] = ...
    learning_errors(G,E1,perm_matrix);
fprintf("accuracy:\nextra edges: %d\nmissing edges: %d\nprecision: " + ...
    "%0.2f\nrecall: %0.2f\nF1 score: %0.2f\n\n\n",extra_edges, ...
    missing_edges, precision, recall, skeleton_F1_score);
       


% 2) finite sample setting. No need to feed the ground truth dag.

MARVEL_Opt2 = {'data', D};
fprintf('finite sample setting with %d samples...\n',num_of_samples);
[E2, Mb_tests, MARVEL_tests] = MARVEL(MARVEL_Opt2{:});
fprintf('MARVEL is done. Number of CI tests: %d \n',Mb_tests+MARVEL_tests);
[extra_edges, missing_edges, precision, recall, skeleton_F1_score] = ....
    learning_errors(G,E2,perm_matrix);
fprintf("accuracy:\nextra edges: %d\nmissing edges: %d\nprecision: " + ...
    "%0.2f\nrecall: %0.2f\nF1 score: %0.2f\n\n\n",extra_edges, ...
    missing_edges, precision, recall, skeleton_F1_score);
