function CI = CI_Test(X,Y,S,D,alpha)

% The CI test function that works for linear Gaussian models
% It checks whether X and Y are independent conditioning on S using the
% data matrix D and significance level alpha

num_of_samples = size(D,1);
c = norminv(1-alpha/2);
DD = D(:,[X,Y,S]);
R =  corrcoef(DD); %cov(DD);
P = inv(R);
ro = -P(1,2)/sqrt(P(1,1)*P(2,2));
zro = 0.5*log((1+ro)/(1-ro));
if abs(zro)<c/sqrt(num_of_samples-size(S,2)-3)
    CI = true;
else
    CI = false;
end
end
