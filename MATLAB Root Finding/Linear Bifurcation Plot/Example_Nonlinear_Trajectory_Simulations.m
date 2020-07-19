%The time values for which we wish to simulate the trajectories of the
%system.
t_values = linspace(0,50,100);

%The number of integral discretisation points.
N=3;

%The initial data.
y0 = @(t) [sin(t),cos(t)];

%The weights and evaluation points to be used for the legendre polynomial
%approximation of our integral.
[legendre_zeros, legendre_weights] = computeGaussLegendreWeights(0, 1, N);

%Produce a trajectory simulation.
[solx, soly] = compute_trajectory_simulation(double(legendre_zeros), double(legendre_weights),t_values,y0);

%Plot the resulting trajectory.
hold on
plot(solx, sin(solx));
plot(solx, cos(solx));
plot(solx, soly(1,:));
plot(solx, soly(2,:));

legend('sin(t)', 'cos(t)', 'x(t)', 'y(t)');
xlabel('t');


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
    dydt = @(t,y,Z) custom_rule_model(delay_times, weights, t,y,Z);
    
    options = odeset('RelTol',1e-7,'Events',@terminalEventFcn);
    sol = dde23(dydt, delay_times, y0, t_values, options);
    
    x = sol.x;
    y = sol.y;
end


%This function uses Gauss-Legendre quadrature to approximate the
%derivative in the system, for a given set of parameters.
function derivative = custom_rule_model(delay_values, weights, t, y, Z)
    integrand = cos(delay_values) .* Z(1,:)';
    integral_approximation = dot(weights, integrand);
    
    derivative_x = y(2);
    derivative_y = -y(1) - y(1) * (sin(2)/4 +1/2) - y(2)*(cos(2)-1)/4 + integral_approximation;
    derivative = [derivative_x, derivative_y]';
end
