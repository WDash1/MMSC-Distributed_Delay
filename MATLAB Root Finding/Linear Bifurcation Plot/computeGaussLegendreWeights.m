%This function may be used to calculate the zeros and integral weights
%associated with the Gauss-Legendre quadrature approximation of an integral
%over a given interval.
function [zeros, weights] = computeGaussLegendreWeights(a, b, n)
    %Compute the required legendre polynomial.
    legendre_polynomial = @(x) legendreP(n,x);
    
    %Compute the corresponding derivative of the polynomial.
    syms x;
    legendre_polynomial_derivative = diff(legendre_polynomial, x, 1);

    %Find the zeros of the polynomial.
    syms y;
    raw_zeros = vpasolve(legendre_polynomial(y) == 0);
    zeros = (raw_zeros.*((b-a)/2)) + (a+b)/2;

    
    %Compute the corresponding weights for the integral we wish to
    %approximate.
    raw_weights = double(2./((1-raw_zeros.^2).*(subs(legendre_polynomial_derivative, raw_zeros).^2)));
    weights = double((b-a)/2 .* raw_weights);
end