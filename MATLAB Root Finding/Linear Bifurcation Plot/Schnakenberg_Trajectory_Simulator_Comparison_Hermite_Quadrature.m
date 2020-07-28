%The time mesh on which we evaluate trajectories.
t_mesh = linspace(100,1200,50000);

%Parameter values
a = 0.5;
b = 5;
tau = 0.05;

%The fixed point for the parameter values.
u_fixed = b + a;
v_fixed = b/((b+a).^2);

%The initial data.
u0 = @(t) u_fixed + 0.001*(sin(sqrt(2)*t) + cos(t));
v0 = @(t) v_fixed + 0.001*(sin(sqrt(2)*t) + cos(t));


%Plot the no delay trajectory.
figure
hold on;
xlabel('{\it u}');
ylabel('{\it v}');
plot(u_fixed,v_fixed,'r*', 'DisplayName','Fixed Point')
[u_values, v_values] = computeNoDelaySchnakenbergTrajectory(a,b,t_mesh,u0(0),v0(0),10e-18,10e9);
plot(u_values,v_values, 'DisplayName',"No delay");

%Plot the fixed delay trajectory.
sol = computeFixedDelaySchnakenbergTrajectory(a,b,tau,max(t_mesh), u0,v0,10e-18,10e9);
%current_y_values = deval(sol, t_mesh);
current_y_values = deval(sol,t_mesh);
u_values = current_y_values(1,:);
v_values = current_y_values(2,:);
plot(u_values,v_values, 'DisplayName',"Fixed delay");
    



%The number of points to use in the quadrature rule.
N=3;

%The sigma values we wish to simulate.
sigma_amt = 0;
sigma_values = linspace(0.001,0.002,sigma_amt);

%Plot the distributed delay trajectories.
for index=1:sigma_amt
    
    
    
    sigma = sigma_values(index);
    [integral_evaluation_point, integral_weights] = computeGaussHermiteWeights(tau, sigma, N);
    
    filtered_evaluation_points = integral_evaluation_point;
    
    sol = computeDistributedDelaySchnakenbergTrajectory2(double(integral_evaluation_point), ...
                                                        double(integral_weights), ...
                                                        a,b,...
                                                        max(t_mesh),...
                                                        u0,v0,...
                                                        10e-9,10e9);
    y_values = deval(sol, t_mesh);
    %y_values = sol.y;
    u_values = y_values(1,:);
    v_values = y_values(2,:);
    
    
    plot(t_mesh,v_values, 'DisplayName',"\sigma = "+string(round(sigma,10)));
    
end
legend
%print('Trajcetory Comparison Plot', '-dpng', '-r300');



function solution = computeDistributedDelaySchnakenbergTrajectory2(...
                                                delay_times, ...
                                                weights, ...
                                                a, ...
                                                b, ...
                                                t_max, ...
                                                u0, ...
                                                v0, ...
                                                relative_tolerance, ...
                                                termination_threshold)

    
    %Generate the derivative function for our system, using the given
    %parameter values.
    dydt = @(t,y,Z) distributedDelaySchnakenberg(weights(:)', a,b,t,y,Z);
    
    %Generate an event function to stop the simulation, if it exceeds the
    %specified threshold.
    eventFunction = @(t,y,Z) terminalEventFcn(termination_threshold,t,y,Z);
    options = odeset('RelTol',relative_tolerance,'Events',eventFunction);

    
    %If zero is an evaluation point in our delay_times vector, then exclude
    %it from the dde23 function call and re-insert it in the underlying
    %distributedDelaySchnakenberg function.
    y0 = @(t) [u0(t),v0(t)];
    if(delay_times(1) == 0)  
        solution = dde23(dydt, delay_times(2:end), y0, [0,t_max], options);
    else
        solution = dde23(dydt, delay_times, y0, [0,t_max], options);
    end
end

function derivative = distributedDelaySchnakenberg(weights, a, b, t,y,Z)
    %If zero was one of the delay values, then include the value of u and v
    %at time t in the vector of function values.
    function_values = Z;
    if(size(Z,2)+1 == size(weights,2))
        function_values = [y, Z];
    end
    
    %Use the provided weights to approximate the integral in the governing
    %equation.
    integrand = ((function_values(1,:)).^2).* function_values(2,:);
    integral_approximation = dot(weights, integrand);

    
    %Return the u and v derivatives.
    derivative_u = a - y(1) + integral_approximation;
    derivative_v = b - integral_approximation;
    derivative = [derivative_u, derivative_v]';
end


%This event function is used to terminate the simulation, if it exceeds a
%suitably large threshold.
function [position,isterminal,direction] = terminalEventFcn(threshold, t,y, Z)
    position = max(abs(y(1)))<threshold; %The value that we want to be zero
    isterminal = 1;  % Halt integration 
    direction = 0; 
end

