%This file is used to approximate the integral of a Gaussian
%distribution, using Gauss-Hermite quadrature.

%The mean and standard deviation of the distribution.
mu = 0.6;
sigma = mu*0.1;

N=31;
[evaluation_points, integrand_weights] = computeGaussHermiteWeights(mu, sigma, N);
%Compute the truncated Gaussian distribution function.
phi = @(x) exp(-(x.^2) ./ 2).* 1/(sqrt(2*pi));
Phi = @(x) (1/2).*(erf(x./sqrt(2)));
K = @(s) 1/sigma .* phi((s - mu) ./ sigma) ./ (Phi(mu / sigma) - Phi(-mu/sigma));
%K = @(s) 1/sigma .* phi((s - mu) ./ sigma);


step_function = @(s) double(s>0) .* double(s<2*mu);

K = @(s) 1/sigma .* phi((s - mu) ./ sigma) .* step_function(s);

truncation_correction_factor = 1/ (Phi(mu / sigma) - Phi(-mu/sigma));
%truncation_correction_factor = 1;

%Plot the truncated Gaussian function and quadrature points.
figure
hold on

s_values = 2.*linspace(0,2*mu, 5000);
plot(s_values, K(s_values));
plot(evaluation_points, K(evaluation_points), '*');


hold off
figure

n_amt=50;
approximation_error = zeros(1,n_amt);
n_values = 1.* (1:n_amt);
for n_index = 1:n_amt
    N=n_values(n_index);
    %The number of quadrature points.
    %N=10;
    [evaluation_points, integrand_weights] = computeGaussHermiteWeights(mu, sigma, N);

    %Approximate the integral of the Gaussian distribution function.
    k_values = K(evaluation_points);
    k_approximation = sum(integrand_weights .* step_function(evaluation_points).*truncation_correction_factor);
    approximation_error(n_index) = abs(k_approximation - 1);
    
    %disp("Approximation to integral of truncated Gaussian distribution (should be 1) ");
    %disp("for N = "+string(N)+" ");
    %disp("approximation is: "+string(k_approximation));
end

hold off
plot(n_values, approximation_error);
xlabel("N");
