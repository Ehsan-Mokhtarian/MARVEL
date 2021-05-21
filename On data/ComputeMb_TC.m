function [Mb,tests] = ComputeMb_TC(D,alpha)
% the "Total Conditioning" algorithm for computing Markov boundaries
% see Pellet and Elisseeff paper

tests = 0;
n = size(D,2);
Mb = zeros(n);

for X=1:(n-1)
    for Y=(X+1):n
        S = [1:(X-1),(X+1):(Y-1),(Y+1):n];
        CI = CI_Test(X,Y,S,D,alpha);
        tests = tests+1;
        if ~CI
            Mb(X,Y)=1;
            Mb(Y,X)=1;
        end
    end
end
end
