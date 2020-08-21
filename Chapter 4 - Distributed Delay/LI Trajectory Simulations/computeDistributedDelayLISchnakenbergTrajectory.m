%@brief                         This function may be used to approximate a
%                               trajectory in the fixed delay Schnakenberg
%                               system, for a given set of parameters.
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
%@param t_max                   This parameter must be a double, which is
%                               >0, that dictates how far in time we wish
%                               to run trajectory simlations for.
%@param u0                      This must be a function which takes a time
%                               value in the interval [-tau, 0] and returns
%                               a value of type double. This dictates the
%                               initial data for the u variable that we
%                               wish to use in our trajectory simulation.
%@param v0                      This must be a function which takes a time
%                               value in the interval [-tau, 0] and returns
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
function solution = computeDistributedDelayLISchnakenbergTrajectory(...
                                                a, ...
                                                b, ...
                                                delay_values, ...
                                                weights, ...
                                                t_max, ...
                                                u0, ...
                                                v0, ...
                                                relative_tolerance, ...
                                                termination_threshold)

                                            
    %Generate an event function to stop the simulation, if it exceeds the
    %specified threshold.
    eventFunction = @(t,y,Z) terminalEventFcn(termination_threshold,t,y,Z);
    options = odeset('RelTol',relative_tolerance,'Events',eventFunction);
    
    
    %Generate the derivative function for our system, using the given
    %parameter values.
    dydt = @(t,y,Z) distributedDelaySchnakenberg(a,b,weights,t,y,Z);
    
    y0 = @(t) [u0(t),v0(t)];
    %Simulate the fixed delay system.
    if(delay_values(1) == 0)
       solution = dde23(dydt, delay_values(2:end), y0, [0,t_max],options);
    else
       solution = dde23(dydt, delay_values, y0, [0,t_max],options);
    end
end


%@brief     This function may be used to compute the derivative of the
%           fixed delay Schnakenberg system, for a given set of parameters
%           and u,v values.
%@param a   This must be a double which is >0.
%@param b   This must be a double which is >0.
%@param t   This must be a double which is >0.
%@param y   This must be a column vector of size 2, where the first
%           entry contains the value of u at time t and the
%           second contains the value of v at time t in the simulation.
%@param Z   This must be a column vector of size 2, where the first entry
%           dictates the value of u at time t-\tau and the second
%           dictates the value of v at time t-\tau.
%@return    This function returns a column vector of size 2, whose first
%           entry dictates the value of the derivative of u at time t and
%           the second dictates the value of the derivative of v at time t.
function derivative = distributedDelaySchnakenberg(a,b,weights, t,y,Z)
    function_values = Z;
    if(size(Z,2)+1 == size(weights,2))
        function_values = [y, Z];
    end
    integrand = ((function_values(1,:)).^2).* function_values(2,:);
    integral_approximation = dot(weights, integrand);
    
    derivative_u = a - y(1) + -2.*y(1).^2.*y(2) +3.*integral_approximation;
    derivative_v = b - y(1).^2.*y(2);
    derivative = [derivative_u, derivative_v]';
    
    %derivative_u = a - y(1) + delay_componnent;
    %derivative_v = b - delay_componnent;
   % derivative_u = a - y(1) + (-2.*(y(1).^2) .* y(2)) +3.*delay_componnent;
   % derivative_v = b - ((y(1).^2) .* y(2));
   % derivative = [derivative_u, derivative_v]';
end


%This event function is used to terminate the simulation, if it exceeds a
%suitably large threshold.
function [position,isterminal,direction] = terminalEventFcn(threshold, t,y, Z)
    position = max(abs(y(1)))<threshold; %The value that we want to be zero
    isterminal = 1;  % Halt integration 
    direction = 0; 
end
