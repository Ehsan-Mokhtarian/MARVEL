clear
clc
addpath('structures/')

%% Load graph and generate data
% See more graphs in 'structures' folder
G = load('mildew').A;
p = size(G,1);
number_of_samples = 50*p;
D = Generate_linear_Gaussian_Data(G, number_of_samples);

%% Calling MARVEL function
alpha = 2/p^2;
Mb = ComputeMb_TC(D, alpha); % Learning  Markov boundaries using TC algorithm
[G_MARVEL, tests, SC] = MARVEL(D, Mb, alpha);
fprintf('MARVEL:\n#CI tests: %d\n',tests);

%% Accuracy of the learned graph
[extra_edges,missing_edges,precision,recall,skeleton_F1_score]=...
    learning_errors(G, G_MARVEL);
fprintf('extra edges: %d\nmissing edges: %d\nprecision: %0.2f\nrecall: %0.2f\nF1 score: %0.2f\n\n\n',extra_edges, missing_edges, precision, recall, skeleton_F1_score);