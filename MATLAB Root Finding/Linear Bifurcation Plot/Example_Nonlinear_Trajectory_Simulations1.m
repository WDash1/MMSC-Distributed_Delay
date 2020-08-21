%The time values for which we wish to simulate the trajectories of the
%system.
t_values = linspace(0,30,100);

%The number of integral discretisation points.
%N=10;

%The initial data.
y0 = @(t) [sin(t),cos(t)];


hold on
%plot(t_values, sin(t_values));
%plot(t_values, cos(t_values));




N=40;

%The weights and evaluation points to be used for the legendre polynomial
%approximation of our integral.
[legendre_zeros, legendre_weights] = computeCompositeSimpson38RuleWeights(0, 1, 1+N*3);

%Produce a trajectory simulation.
[solx, soly] = compute_trajectory_simulation(double(legendre_zeros), double(legendre_weights),t_values,y0);

%cauchy_sequence(N) = max(abs(soly(1,:) - sin(solx)));
cauchy_sequence(N) = max(abs(soly(2,:) - cos(solx)));
    
%Plot the resulting trajectory.
plot(solx, soly(1,:), 'DisplayName','x(t)');
plot(solx, soly(2,:), 'DisplayName','y(t)');

xlabel('t');
legend



%This event function is used to terminate the simulation, if it exceeds a
%suitably large threshold.
function [position,isterminal,direction] = terminalEventFcn(t,y, Z)
    position = max(abs(y(1)))<10^5; % The value that we want to be zero
    isterminal = 1;  % Halt integration 
    direction = 0; 
end


%This function produces a trajectory simulation, for a given set of
%parameters.
function [x,y] = compute_trajectory_simulation(delay_times, weights, t_values, y0)
    dydt = @(t,y,Z) custom_rule_model(delay_times(:)', weights(:)', t,y,Z);
    
    options = odeset('RelTol',1e-10,'Events',@terminalEventFcn);
    
    if(delay_times(1) == 0)  
        sol = dde23(dydt, delay_times(2:end), y0, t_values, options);
    else
        sol = dde23(dydt, delay_times, y0, t_values, options);
    end
    
    x = sol.x;
    y = sol.y;
end


%This function uses Gauss-Legendre quadrature to approximate the
%derivative in the system, for a given set of parameters.
function derivative = custom_rule_model(delay_values, weights, t, y, Z)



    function_values = Z;
    if(size(Z,2)+1 == size(delay_values,2))
        function_values = [y, Z];
    end
    
    integrand = cos(delay_values) .* sqrt(1-((function_values(1,:)).^2)).* sign(cos(t-delay_values));
    integral_approximation = dot(weights, integrand);
    
    derivative_x = y(2);
    derivative_y = -y(1) - y(2) * (sin(2)/4 +1/2) - y(1)*(1-cos(2))/4 + integral_approximation;
    derivative = [derivative_x, derivative_y]';
end



