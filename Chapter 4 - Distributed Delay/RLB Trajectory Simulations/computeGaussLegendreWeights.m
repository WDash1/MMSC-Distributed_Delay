%@brief     This function may be used to calculate the zeros and integral
%           weights associated with the Gauss-Legendre quadrature
%           approximation of an integral over a given interval.
%@param a   This parameter dictates the lower limit of the integral
%           that we wish to use the Gauss-Legendre scheme to approximate.
%           This argument must be a double, which is <b.
%@param b   This parameter dictates the upper limit of the integral
%           that we wish to use the Gauss-Legendre scheme to approximate.
%           This argument must be a double, which is >a.
%@param n   This argument dictates the number of evaluation points
%           that we wish to use to generate for approximating our
%           integral. This argument must be an odd integer which is >0.
%@return    This function returns the evaluation points to be used for
%           approximating our desired integral, as well as the
%           corresponding weight for each of these evaluation points.
function [evaluation_points, weights] = computeGaussLegendreWeights(a, b, n)


    %Compute the required legendre polynomial.
    legendre_polynomial = @(x) legendreP(n,x);
    
    %Compute the corresponding derivative of the polynomial.
    syms x;
    legendre_polynomial_derivative = diff(legendre_polynomial, x, 1);

    %Find the zeros of the polynomial and scale them to lie in our
    %integration domain.
    syms y;
    raw_zeros = vpasolve(legendre_polynomial(y) == 0);
    zeros = raw_zeros.*((b-a)/2) + (a+b)/2;
    evaluation_points = double(zeros);
    
    %Compute the corresponding weights for the integral we wish to
    %approximate and scale them to our integral domain.
    raw_weights = (2./((1-raw_zeros.^2).*(subs(legendre_polynomial_derivative, raw_zeros).^2)));
    weights = double((b-a)/2 .* raw_weights);
end