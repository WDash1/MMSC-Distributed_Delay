%This file is used to approximate the integral of a Gaussian
%distribution, using Gauss-Hermite quadrature.

%The mean and standard deviation of the distribution.
mu = 0.06;
sigma = 0.03;

%The number of quadrature points.
N=3;
[evaluation_points, integrand_weights] = computeGaussHermiteWeights(mu, sigma, N);

%Compute the truncated Gaussian distribution function.
phi = @(x) exp(-(x.^2) ./ 2).* 1/(sqrt(2*pi));
Phi = @(x) (1/2).*(erf(x./sqrt(2)));
K = @(s) 1/sigma .* phi((s - mu) ./ sigma) ./ (Phi(mu / sigma) - Phi(-mu/sigma));
%K = @(s) 1/sigma .* phi((s - mu) ./ sigma);

%Plot the truncated Gaussian function and quadrature points.
figure
hold on
s_values = 2*mu.*linspace(0,1, 5000);
plot(s_values, K(s_values));
plot(evaluation_points, K(evaluation_points), '*');


%Approximate the integral of the Gaussian distribution function.
k_values = K(evaluation_points);
k_approximation = sum(integrand_weights);
disp("Approximation to integral of truncated Gaussian distribution (should be 1) ");
disp("for N = "+string(N)+" ");
disp("approximation is: "+string(k_approximation));
