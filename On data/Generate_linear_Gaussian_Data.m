function D = Generate_linear_Gaussian_Data(G, number_of_samples)
% Variances are random and between 1 to 3
% Coefficients are random and in [-1,0.5] \cup [0.5,1]

p = size(G,1);
N = randn(number_of_samples,p)*diag(1 + 0.732*rand(1,p));
A = (0.5+0.5*rand(p)) .* ((-1).^(rand(p)>0.5));
A = G.*A;
D = N/(eye(p)-A);
end
