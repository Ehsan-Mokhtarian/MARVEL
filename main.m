close all
clear
clc
addpath('./structures');

%% Configuration 

num_of_samples = 2000;

% select ground-truth graph
data_structure  =  ["cancer","asia","sachs","insurance","water"...
    ,"mildew","alarm","barley","hailfinder","carpo","hepar2"...
    ,"win95pts","diabetes","pathfinder","andes","pigs","Link"...
    ,"munin2","munin4","munin1","munin3"];
selected_structure = 12;
fprintf('selected structure: %s\n\n', data_structure(selected_structure));

% generate data
[G,D,perm_matrix] = Generate_Data(data_structure(selected_structure), num_of_samples);

%% How to call MARVEL

%  MARVEL Options:

% 'alpha': the parameter of Fischer's Z transform
% 'initial_mb': the algorithm to discover the Markov blanket info.
% algorithms: 'TC', 'GS', 'TCbw'
% 'check_Mb_inclusion': true/false, is preprocess done?
% 'top_down_subsets' true/false: the order in which MARVEL checks the
% subsets of a set, small to large or vice-versa
% 'use_L0': true/false save the CI test results
% 'cond_indep': CI tests, 'oracle' or 'stat' (data)
% 'G': % Ground truth 



% two examples:

% 1) oracle setting. Ground truth, G, must be fed to do the CI tests.

MARVEL_Opt1 = {'check_Mb_inclusion',false,'top_down_subsets', true,...
    'use_L0', true, 'cond_indep', 'oracle', 'G', G};
fprintf('orcale setting...\n');
[E1, mb_tests, marvel_tests] = MARVEL(MARVEL_Opt1{:});
fprintf('MARVEL is done. Number of CI tests: %d \n',mb_tests+marvel_tests);
[extra_edges, missing_edges, precision, recall, skeleton_F1_score] = learning_errors(G,E1,perm_matrix);
fprintf('accuracy:\nextra edges: %d\nmissing edges: %d\nprecision: %0.2f\nrecall: %0.2f\nF1 score: %0.2f\n\n\n',extra_edges, missing_edges, precision, recall, skeleton_F1_score);
       


% 2) finite sample setting. No need to feed the ground truth dag.

MARVEL_Opt2 = {'data', D, 'alpha',0.01, 'initial_mb', 'TC', 'cond_indep', 'stat'};
fprintf('finite sample setting with %d samples...\n',num_of_samples);
[E2, mb_tests, marvel_tests] = MARVEL(MARVEL_Opt2{:});
fprintf('MARVEL is done. Number of CI tests: %d \n',mb_tests+marvel_tests);
[extra_edges, missing_edges, precision, recall, skeleton_F1_score] = learning_errors(G,E2,perm_matrix);
fprintf('accuracy:\nextra edges: %d\nmissing edges: %d\nprecision: %0.2f\nrecall: %0.2f\nF1 score: %0.2f\n\n\n',extra_edges, missing_edges, precision, recall, skeleton_F1_score);


