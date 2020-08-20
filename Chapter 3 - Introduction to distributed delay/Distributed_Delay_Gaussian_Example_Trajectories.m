t_values = linspace(0,1,500);

u0 = @(t) exp(t);
v0 = @(t) exp(t);
y0 = @(t) [u0(t),v0(t)];


mu = 0.1;
sigma = 0.1;


N = 40;    

[gauss_legendre_zeros, gauss_legendre_weights] = computeGaussLegendreWeights(0, 1, N);
[~, legendre_soly] = compute_trajectory_simulation(gauss_legendre_zeros, ...
                                                    gauss_legendre_weights, ...
                                                    0, ...
                                                    2*mu, ...
                                                    mu, ...
                                                    sigma, ...
                                                    t_values, ...
                                                    y0);
    
hold on
plot(t_values, legendre_soly(1,:), '-', 'DisplayName', 'x(t)');
plot(t_values, legendre_soly(2,:), '-', 'DisplayName', 'y(t)');
plot(t_values, exp(t_values), '-', 'DisplayName', 'exp(t)')
xlabel('{\it t}', 'FontSize', 20);
ylabel('{\it y(t)}', 'FontSize', 20);

legend
ax = gca;
ax.FontSize = 20; 
filename = "Distributed_Delay_Example_Gaussian_Trajectory.eps";
print('-depsc', '-tiff', '-r300', '-painters', filename);

function [x,y] = compute_trajectory_simulation(delay_times,weights, a,b,mu,sigma, t_values, y0)
    delays = delay_times;
    if(delay_times(1) == 0)
        delays = delay_times(2:end);
    end
    dydt = @(t,y,Z) linearDerivativeExample(delay_times,weights,a,b,mu,sigma, t,y,Z);
    options = odeset('RelTol',10e-15,'AbsTol',10e-15);
    sol = dde23(dydt, delays, y0, t_values,options);
    
    x = sol.x;
    y = deval(sol, t_values);
    
end



function derivative = linearDerivativeExample(delay_times,weights, a,b,mu,sigma, t,y,Z)
    function_values = Z;
    if(size(Z,2)+1 == size(weights,2))
        function_values = [y, Z];
    end
    
    
    integrand = exp(-((mu-delay_times(:)').^2)/2).*function_values(2,:).^2;
    integral_approximation = dot(weights, integrand);
    
    
    
    mu1 = (2*t*sigma - mu) / (2*sigma.^2 - 1);
    sigma1 = sigma/(sqrt(1-2*sigma^2));
    
    
    derivative_x = y(2);
    derivative_y = y(1) + integral_approximation - ...
                    (exp(t.^2 - mu.^2/(2.*sigma.^2))./(sigma.*sqrt(2*pi))) .* ...
                    exp(-((2*sigma.^2 - 1)/(2.*sigma)).*(((mu - 2*t*sigma.^2)./(2*sigma.^2 - 1)).^2)) .* ...
                    sqrt(2*pi)*sigma1*(normcdf(b, mu1, sigma1)- normcdf(a,mu1,sigma1));
                    %sqrt(2*pi)*sigma1*(normcdf(b, mu1, sigma1)- normcdf(a,mu1,sigma1));
    
    
    
                
     %derivative_y = y(1) + integral_approximation - y(2).^2 .* ...
     %               exp((((2*(sigma^2) - mu).^2)-mu.^2)/(sigma.^2)) .* ...
     %               (normcdf(b+2*sigma.^2, mu, sigma) - normcdf(a+2*sigma.^2, mu, sigma));
    
     derivative = [double(derivative_x), double(derivative_y)]';
end