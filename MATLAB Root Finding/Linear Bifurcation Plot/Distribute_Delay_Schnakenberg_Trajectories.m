%The time values for which we wish to simulate the trajectories of the
%system.
t_values = linspace(0,100,100);

%The number of integral discretisation points.
%N=10;

%The initial data.
y0 = @(t) [sin(t),cos(t)];


hold off
%plot(t_values, sin(t_values));
%plot(t_values, cos(t_values));


param_amt = 4;
a_values = [0.01, 2, 3.5, 5];
b_values = [0.5, 2, 3, 5];
tau = 0.03;
a = 0.01;
b = 0.5;
sigma = 0.001;

n_amt = 10;
n_values = ((1:n_amt).*2);

hold on;
parfor index = 1:param_amt
    current_a = a_values(index);
    current_b = b_values(index);

    cauchy_sequence = generateCauchySequence(current_a,current_b,tau,sigma,n_values, 5);
    
    
    figure
    plot(n_values(2:end), cauchy_sequence);
    title("a="+string(current_a)+" b="+string(current_b));
    xlabel('N');
    print('Cauchy_Plot'+string(index), '-dpng', '-r300');

end




function cauchy_sequence = generateCauchySequence(a,b,tau,sigma,n_values, t_max)
    cauchy_sequence = zeros(1,size(n_values,2)-1);
    
    
    comparison_mesh = linspace(0, t_max, 100);
    
    u_fixed = b+a;
    v_fixed = b/((b+a).^2);

    
    u0 = @(t) u_fixed + 1*(sin(sqrt(2)*t) + cos(t));
    v0 = @(t) v_fixed + 1*(sin(sqrt(2)*t) + cos(t));
    y0 = @(t) [u0(t), v0(t)];


    last_u_values = 0;
    last_v_values = 0;
    
    for index = 1:size(n_values,2)
        current_n = n_values(index);
        [legendre_zeros, legendre_weights] = computeGaussLegendreWeights(0, 2*tau, current_n);
        sol = compute_trajectory_simulation(a,b,sigma,double(legendre_zeros), double(legendre_weights),t_max,y0);
        

        new_comparison_mesh = comparison_mesh;
        if(max(sol.x)<max(comparison_mesh))
            new_comparison_mesh = comparison_mesh .* (comparison_mesh < max(sol.x));
        end
        current_y_values = deval(sol, new_comparison_mesh);
        current_u_values = current_y_values(1,:);
        current_v_values = current_y_values(2,:);
        
        if(index>1)
            max_u_difference = max(abs(current_u_values - last_u_values));
            max_v_difference = max(abs(current_v_values - last_v_values));
            cauchy_sequence(index-1) = max(max_u_difference, max_v_difference);
        end
        
        last_u_values = current_u_values;
        last_v_values = current_v_values;
    end
end


%This event function is used to terminate the simulation, if it exceeds a
%suitably large threshold.
function [position,isterminal,direction] = terminalEventFcn(t,y, Z)
    position = max(abs(y(1)))<10^9; % The value that we want to be zero
    isterminal = 1;  % Halt integration 
    direction = 0; 
end


%This function produces a trajectory simulation, for a given set of
%parameters.
function solution = compute_trajectory_simulation(a,b,sigma,delay_times, weights, t_max, y0)

    tau = max(delay_times)/2;

    phi = @(x) exp(-(x.^2) ./ 2).* 1/(sqrt(2*pi));
    Phi = @(x) (1/2).*(erf(x./sqrt(2)));

    
    K = @(s) 1/sigma .* phi((s - tau) ./ sigma) ./ (Phi(tau / sigma) - Phi(-tau/sigma));
    
    kernel_values = K(delay_times);
    
    dydt = @(t,y,Z) custom_rule_model(delay_times(:)', weights(:)', a,b,kernel_values(:)',t,y,Z);
    
    options = odeset('RelTol',1e-10,'Events',@terminalEventFcn);
    
    if(delay_times(1) == 0)  
        sol = dde23(dydt, delay_times(2:end), y0, [0,t_max], options);
    else
        sol = dde23(dydt, delay_times, y0, [0,t_max], options);
    end
    
    solution = sol;
end


%This function uses Gauss-Legendre quadrature to approximate the
%derivative in the system, for a given set of parameters.
function derivative = custom_rule_model(delay_values, weights, a, b, kernel_values,t, y, Z)


    function_values = Z;
    if(size(Z,2)+1 == size(delay_values,2))
        function_values = [y, Z];
    end
    
    integrand = ((function_values(1,:)).^2).* function_values(2,:) .* kernel_values;
    integral_approximation = dot(weights, integrand);
    
    derivative_u = a - y(1) + integral_approximation;
    derivative_v = b - integral_approximation;
    derivative = [derivative_u, derivative_v]';
end



