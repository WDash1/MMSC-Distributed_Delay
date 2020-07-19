%The time values for which we wish to simulate the trajectories of the
%system.
t_values = linspace(0,50,100);

%The number of integral discretisation points.
N=3;

%The parameter values that we wish to use for our simulation.
alpha_amt = 10;
beta_amt = 10;
%alpha_values = linspace(-2.5,-1.75,alpha_amt);
%beta_values = linspace(3.25,3.75,beta_amt);
alpha_values = linspace(-10,2,alpha_amt);
beta_values = linspace(-20,10,beta_amt);

%The initial data.
y0 = @(t) [sin(t),cos(t)];


[legendre_zeros, legendre_weights] = computeGaussLegendreWeights(0, 1, N);


[solx, soly] = compute_trajectory_simulation(double(legendre_zeros), double(legendre_weights), N,t_values,y0);

hold on
plot(solx, soly(1,:));
plot(solx, soly(2,:));



legend('x', 'y');
xlabel('t');

    



function [position,isterminal,direction] = terminalEventFcn(t,y, Z)
    position = max(abs(y(1)))<10^3; % The value that we want to be zero
    isterminal = 1;  % Halt integration 
    direction = 0; 
end




%This function produces a trajectory simulation, for a given set of
%parameters.
function [x,y] = compute_trajectory_simulation(delay_times, weights, N, t_values, y0)
%    delay_times = linspace(0,1,N);
%    delay_times = delay_times(2:end);
    
    %dydt = @(t,y,Z) simpsons_38_rule_model( N, t,y,Z);
    dydt = @(t,y,Z) custom_rule_model(weights, t,y,Z);
    
    options = odeset('RelTol',1e-7,'Events',@terminalEventFcn);
    sol = dde23(dydt, delay_times, y0, t_values, options);
    
    x = sol.x;
    y = sol.y;
end



%This function uses the trapezium rule to approximate the derivative in the
%system, for a given set of parameters.
function derivative = trapezium_rule_model( N, t,y,Z)
    sum_values = sum(Z(1,:)) - ((y(1)+Z(1,end))/2);
    derivative_x = y(2);
    derivative_y = -y(1) + 1/N * sum_values - sin(1)*y(1) + (1-cos(1))*y(2);
    derivative = [derivative_x, derivative_y]';
end

%This function uses the composite simpsons rule to approximate the
%derivative in the system, for a given set of parameters.
function derivative = simpsons_rule_model( N, t,y,Z)
    even_indicies = linspace(3, N-2, ((N-1)/2)-1);
    odd_indicies = linspace(2, N-1, ((N-1)/2));
    
    sum_value = y(1)+2*sum(Z(1,even_indicies)) + 4*sum(Z(1,odd_indicies))+Z(1,end);
    derivative_x = y(2);
    derivative_y = -y(1) + 1/(3*N) * sum_value - sin(1)*y(1) + (1-cos(1))*y(2);
    derivative = [derivative_x, derivative_y]';
end

function derivative = custom_rule_model(weights, t, y, Z)
    integrand = Z(1,:);
    integral_approximation = dot(weights, integrand);
    
    derivative_x = y(2);
    derivative_y = -y(1) - sin(1)*y(1) + (1-cos(1))*y(2) + integral_approximation;
    derivative = [derivative_x, derivative_y]';
end

%This function uses the composite simpsons 38 rule to approximate the
%derivative in the system, for a given set of parameters.
function derivative = simpsons_38_rule_model( N, t,y,Z)
    third_indicies = linspace(4,N-3, ((N-1)/3)-1);
    
    first_third_indicies = linspace(3,N-1, ((N-1)/3));
    second_third_indicies = linspace(3,N-1, ((N-1)/3));
    
    
    
    sum_value = y(1) + 3*(sum(Z(1,first_third_indicies))+sum(Z(1,second_third_indicies)))...
                +2*sum(Z(1,third_indicies)) + Z(1,end);
            
            
            
    derivative_x = y(2);
    derivative_y = -y(1) + 3/(8*N) * sum_value - sin(1)*y(1) + (1-cos(1))*y(2);
    derivative = [derivative_x, derivative_y]';
end