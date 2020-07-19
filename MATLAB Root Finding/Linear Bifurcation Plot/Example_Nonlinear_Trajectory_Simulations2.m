%The time values for which we wish to simulate the trajectories of the
%system.
t_values = linspace(0,3,100);

%The number of integral discretisation points.
N=3;

%The initial data.
y0 = @(t) exp(t);

%The weights and evaluation points to be used for the legendre polynomial
%approximation of our integral.
[legendre_zeros, legendre_weights] = computeGaussLegendreWeights(0, 1, N);

%Produce a trajectory simulation.
[solx, soly] = compute_trajectory_simulation(double(legendre_zeros), double(legendre_weights), t_values,y0);

%Plot the resulting trajectory.
hold on
plot(solx, exp(solx));
plot(solx, soly);

legend('exp(t)', 'x(t)');

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


%This function computes the derivative at each step in the simulation, for
%our system.
function derivative = custom_rule_model(delay_values, weights, t, y, Z)
    integrand = exp(5.*delay_values) .* ((Z').^4);
    integral_approximation = dot(weights, integrand);
    
    derivative = y + (y.^4)*(1-exp(1)) + integral_approximation;
end

