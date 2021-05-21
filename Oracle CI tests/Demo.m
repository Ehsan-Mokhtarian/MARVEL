clear
clc
addpath('structures/')

%% Load Graph
% See more graphs in 'structures' folder
G = load('diabetes').A;

%% Calling MARVEL function

Mb = ComputeMb_oracle(G); % The true Markov boundaries
[G_MARVEL, tests, SC] = MARVEL_oracle(G, Mb);
