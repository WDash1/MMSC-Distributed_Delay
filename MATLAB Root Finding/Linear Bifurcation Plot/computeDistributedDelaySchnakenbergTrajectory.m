%@brief                         This function may be used to approximate a
%                               trajectory in the distributed delay
%                               Schnakenberg system, for a given set of
%                               parameters and integrand evaluation points
%                               and weights.
%@param delay_times             This must be a column vector of type
%                               double, which is the same size as weights
%                               and dictates the delay time values that we
%                               wish to use to approximate the integrand
%                               term in the governing equation.
%@param weights                 This must be a column vector of type
%                               double, which is the same size as
%                               delay_times and dictates the weights of
%                               each of the integrand evaluation points in
%                               the delay_times vector, that we wish to use
%                               to approximate the integral term in the
%                               governing equation.
%@param a                       This parameter must be a double, which is
%                               >0, that dictates the value of the a
%                               parameter to be used in simulations of the
%                               model.
%@param b                       This parameter must be a double, which is
%                               >0, that dictates the value of the b
%                               parameter to be used in simulations of the
%                               model.
%@param tau                     This parameter must be a double, which is
%                               >0, that dictates the time delay value
%                               that we wish to use in simulations of
%                               the model.
%@param sigma                   This parameter must be a double, which is
%                               >0, that dictates the standard deviation of
%                               the truncated normal distribution that we
%                               wish to use as the kernel function for
%                               trajectory simulations.
%@param t_max                   This parameter must be a double, which is
%                               >0, that dictates how far in time we wish
%                               to run trajectory simlations for.
%@param u0                      This must be a function which takes a time
%                               value in the interval 
%                               [-max(delay_times), 0] and returns
%                               a value of type double. This dictates the
%                               initial data for the u variable that we
%                               wish to use in our trajectory simulation.
%@param v0                      This must be a function which takes a time
%                               value in the interval 
%                               [-max(delay_times), 0] and returns
%                               a value of type double. This dictates the
%                               initial data for the v variable that we
%                               wish to use in our trajectory simulation.
%@param relative_tolerance      This must be a value of type double, that
%                               we wish to use as the relative tolernace
%                               for the trajectory simulation.
%@param termination_threshold   This must a value of type double, which
%                               dictates the maximum absolute value that
%                               may be obtained by either u or v, before
%                               the simulation is halted.
%@return                        This function returns a struct with fields:
%                               x,y,yp and solver. The values of u and v 
%                               may be evaluated on a desired mesh by using
%                               the deval function.
function solution = computeDistributedDelaySchnakenbergTrajectory(...
                                                delay_times, ...
                                                weights, ...
                                                a, ...
                                                b, ...
                                                sigma, ...
                                                t_max, ...
                                                u0, ...
                                                v0, ...
                                                relative_tolerance, ...
                                                termination_threshold)

    
  
  
    
    %Produce a symmetric truncated Gaussian distribution, for the kernl
    %function.
    tau = max(delay_times)/2;
    phi = @(x) exp(-(x.^2) ./ 2).* 1/(sqrt(2*pi));
    Phi = @(x) (1/2).*(erf(x./sqrt(2)));
    K = @(s) 1/sigma .* phi((s - tau) ./ sigma) ./ (Phi(tau / sigma) - Phi(-tau/sigma));
    kernel_values = K(delay_times);

    
    %Generate the derivative function for our system, using the given
    %parameter values.
    dydt = @(t,y,Z) distributedDelaySchnakenberg(weights(:)', a,b,kernel_values(:)',t,y,Z);
    
    %Generate an event function to stop the simulation, if it exceeds the
    %specified threshold.
    eventFunction = @(t,y,Z) terminalEventFcn(termination_threshold,t,y,Z);
    options = odeset('RelTol',relative_tolerance,'AbsTol',relative_tolerance,'Events',eventFunction);

    
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


%@brief                 This function may be used to compute the derivative
%                       of the distributed delay Schnakenberg system, for a
%                       given set of parameters, integrand weights and u,v
%                       values.
%@param weights         This must be a column vector of type double, that
%                       dictates the weighting of each of the fixed delay 
%                       terms that we wish to use to our approximation of 
%                       the integral term in the distributed delay equation.
%@param a               This must be a double which is >0.
%@param b               This must be a double which is >0.
%@param kernel_values   This must be a column vector of type double, that
%                       dictates the value of the kernel function at each
%                       of the function evaluation points in the integral
%                       approximation.
%@param t               This must be a double which is >0.
%@param y               This must be a column vector of size 2, where the
%                       first entry contains the value of u at time t and
%                       the second contains the value of v at time t in
%                       the simulation.
%@param Z               This must be a matrix of size column vector of size
%                       2xn, where the (1,i)th value corresponds
%                       to the value of u at time t-\tau_{i}, where 
%                       \tau_{i} is used to denote the ith delay value used
%                       to approximate the integral term in the governing
%                       equation.
%@return                This function returns a column vector of size 2,
%                       whose first entry dictates the value of the
%                       derivative of u at time t and the second dictates
%                       the value of the derivative of v at time t.
function derivative = distributedDelaySchnakenberg(weights, a, b, kernel_values,t,y,Z)
    %If zero was one of the delay values, then include the value of u and v
    %at time t in the vector of function values.
    function_values = Z;
    if(size(Z,2)+1 == size(kernel_values,2))
        function_values = [y, Z];
    end
    
    %Use the provided weights to approximate the integral in the governing
    %equation.
    integrand = ((function_values(1,:)).^2).* function_values(2,:) .* kernel_values;
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
