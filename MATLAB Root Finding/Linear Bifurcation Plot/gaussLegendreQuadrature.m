function integral_approximation = gaussLegendreQuadrature(a, b, f, n)
    %Map the function onto the interval [-1, 1]
    translated_f = @(x) (b-a)/2 * f(((b-a)*x)/2 + (a+b)/2);

    %Compute the required legendre polynomial.
    legendre_polynomial = @(x) legendreP(n,x);
    
    %Compute the corresponding derivative of the polynomial.
    syms x;
    legendre_polynomial_derivative = diff(legendre_polynomial, x, 1);
    
    %Find the zeros of the polynomial.
    syms y;
    zeros = vpasolve(legendre_polynomial(y) == 0);

    %Compute the corresponding weights for the integral we wish to
    %approximate.
    weights = 2./((1-zeros.^2).*(subs(legendre_polynomial_derivative, zeros).^2));
    
    %Evaluate the input function at the zeros of the legendre polynomial.
    function_values = arrayfun(translated_f, zeros);
    
    %Compute the approximation of the integral.
    integral_approximation = dot(function_values, weights);
end