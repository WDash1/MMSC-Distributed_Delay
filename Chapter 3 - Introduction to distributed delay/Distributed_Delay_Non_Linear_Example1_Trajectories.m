t_values = linspace(0,25,500);

y0 = @(t) exp(t);


N = 20;    

[gauss_legendre_zeros, gauss_legendre_weights] =  computeTrapeziumRuleWeights(0, 1, N);
[~, legendre_soly] = compute_trajectory_simulation(gauss_legendre_zeros, ...
                                                    gauss_legendre_weights, ...
                                                    t_values, ...
                                                    y0);
    
hold on
plot(t_values, exp(t_values),'-','DisplayName','exp(t)');
plot(t_values, legendre_soly, '-', 'DisplayName', 'Gauss Legendre');
xlabel('{\it t}', 'FontSize', 20);
ylabel('{\it y(t)}', 'FontSize', 20);

legend
ax = gca;
ax.FontSize = 20; 
filename = "Distributed_Delay_Example1_Non_Linear_Trajectory.eps";
print('-depsc', '-tiff', '-r300', '-painters', filename);

function [x,y] = compute_trajectory_simulation(delay_times,weights, t_values, y0)
    delays = delay_times;
    if(delay_times(1) == 0)
        delays = delay_times(2:end);
    end
    dydt = @(t,y,Z) linearDerivativeExample(delay_times,weights, t,y,Z);
    options = odeset('RelTol',10e-10,'AbsTol',10e-10);
    sol = dde23(dydt, delays, y0, t_values,options);
    
    x = sol.x;
    y = deval(sol, t_values);
    
end



function derivative = linearDerivativeExample(delay_times, weights, t,y,Z)
    function_values = Z;
    if(size(Z,2)+1 == size(weights,2))
        function_values = [y, Z];
    end
    integrand = (exp(5.*delay_times(:)) .* ((function_values').^4));
    integral_approximation = dot(weights(:), integrand);
    derivative = y +(y.^4)*(1-exp(1))+integral_approximation;
end