function [evaluation_points, weights] = computeGaussHermiteWeights(mu, sigma, n)

    %Compute the required hermite polynomial.
    hermite_polynomial = @(x) hermiteH(n,x);
   
    %Find the zeros of the polynomial and scale them to lie in our
    %integration domain.
    syms y;
    raw_zeros = vpasolve(hermite_polynomial(y) == 0);
    zeros = sqrt(2).*sigma.*raw_zeros + mu;
    evaluation_points = double(zeros);
    
    lower_order_hermite = @(x) hermiteH(n-1,x);
    
    raw_weights = (2^(n-1) * factorial(n) * sqrt(pi)) ./ ((n^2) .* (lower_order_hermite(raw_zeros).^2));
    weights = raw_weights ./ (sqrt(pi));
    
end