%The time mesh on which we evaluate trajectories.
t_mesh = linspace(0,10,500);

%Parameter values
a = 0.01;
b = 0.5;


tau = 0.06
%The fixed point for the parameter values.
u_fixed = b + a;
v_fixed = b/((b+a).^2);

%The initial data.
u0 = @(t) u_fixed + 0.001*(sin(sqrt(2)*t) + cos(t));
v0 = @(t) v_fixed + 0.001*(sin(sqrt(2)*t) + cos(t));


%The different sigma values we wish to simulate.
sigma_amt = 5;
sigma_values = linspace(0.01,0.03,sigma_amt);

%The different n values we wish to use for constructing our Cauchy
%sequence.
n_amt = 20;
n_values = 1.*(1:n_amt);

hold on
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')


for sigma_index = 1:sigma_amt
    current_sigma = sigma_values(sigma_index);
    cauchy_sequence = zeros(1,n_amt);

    last_u_values = zeros(size(t_mesh));
    last_v_values = zeros(size(t_mesh));

    current_u_values = last_u_values;
    current_v_values = last_v_values;

    for n_index = 1:n_amt
        current_n = n_values(n_index);
        [integral_evaluation_point, integral_weights] = computeGaussHermiteWeights(tau, current_sigma, current_n);
    
        sol = computeDistributedDelaySchnakenbergTrajectory2(double(integral_evaluation_point), ...
                                                            double(integral_weights), ...
                                                            a,b,...
                                                            max(t_mesh),...
                                                            u0,v0,...
                                                            10e-9,10e9);
        y_values = deval(sol, t_mesh);
        current_u_values = y_values(1,:);
        current_v_values = y_values(2,:);
    
        difference = max(sqrt((last_u_values - current_u_values).^2+(last_v_values - current_v_values).^2));
        cauchy_sequence(n_index) = difference;
    
        last_u_values = current_u_values;
        last_v_values = current_v_values;
    end
    loglog(n_values, cauchy_sequence, 'DisplayName',"\sigma = "+string(round(current_sigma,10)));
    
end


legend
xlabel('N');
%print('Gauss_Hermite_Cauchy_Curve_Plot', '-dpng', '-r300');


%This function computes a simulation of trjectory in the distributed delay
%version of the Schnakenberg equation.
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

