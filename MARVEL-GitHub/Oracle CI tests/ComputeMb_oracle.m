function [Mb, tests] = ComputeMb(G)
% Compute the Markov boundaries in the oracle setting
n = size(G,1);
tests = n*(n-1)/2;
Mb = G+G';
for X=1:n
    Ch = find(G(X,:));
    for Y=Ch
        Pa = G(:,Y);
        Pa(X)=0;
        Mb(X,Pa==1)=1;
    end
end
end