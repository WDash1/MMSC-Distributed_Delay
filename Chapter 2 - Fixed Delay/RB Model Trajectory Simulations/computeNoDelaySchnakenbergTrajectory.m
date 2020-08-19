%@brief                         This function may be used to approximate a
%                               trajectory in the no delay Schnakenberg
%                               system, for a given set of parameters.
%@param a                       This parameter must be a double, which is
%                               >0, that dictates the value of the a
%                               parameter to be used in simulations of the
%                               model.
%@param b                       This parameter must be a double, which is
%                               >0, that dictates the value of the b
%                               parameter to be used in simulations of the
%                               model.
%@param t_mesh                  This parameter dictates the time values at
%                               which we wish to know the value of u and v
%                               in our trajectory simulation.
%@param u0                      This must be a value of type double, that
%                               dictates the initial condition for u that
%                               we wish to use in our simulation.
%@param v0                      This must be a value of type double, that
%                               dictates the initial condition for v that
%                               we wish to use in our simulation.
%@param relative_tolerance      This must be a value of type double, that
%                               we wish to use as the relative tolernace
%                               for the trajectory simulation.
%@param termination_threshold   This must a value of type double, which
%                               dictates the maximum absolute value that
%                               may be obtained by either u or v, before
%                               the simulation is halted.
%@return                        This function returns the column vectors
%                               u_values and v_values. These vectors
%                               respectively correspond to the value of u
%                               and v evaluated on the t_mesh argument
%                               provided.
function [u_values, v_values] = computeNoDelaySchnakenbergTrajectory(...
                                                a, ...
                                                b, ...
                                                t_mesh, ...
                                                u0, ...
                                                v0, ...
                                                relative_tolerance, ...
                                                termination_threshold)

    %Generate the derivative function for our system, using the given
    %parameter values.
    dydt = @(t,y) noDelaySchnakenberg(a,b,t,y);
    
    %Generate an event function to stop the simulation, if it exceeds the
    %specified threshold.
    eventFunction = @(t,y) terminalEventFcn(termination_threshold,t,y);
    options = odeset('RelTol',relative_tolerance,'Events',eventFunction);
    
    %Simulate the trajectory on the mesh and return the result.
    [~, y_values] = ode15s(dydt, t_mesh, [u0,v0], options);
    u_values = y_values(:,1);
    v_values = y_values(:,2);
end


%@brief     This function may be used to compute the derivative of the
%           no delay Schnakenberg system, for a given set of parameters
%           and u,v values.
%@param a   This must be a double which is >0.
%@param b   This must be a double which is >0.
%@param t   This must be a double which is >0.
%@param y   This must be a column vector of size 2, where the first
%           entry contains the value of u at time t and the
%           second contains the value of v at time t in the simulation.
%@return    This function returns a column vector of size 2, whose first
%           entry dictates the value of the derivative of u at time t and
%           the second dictates the value of the derivative of v at time t.
function derivative = noDelaySchnakenberg(a,b,t,y)
    nonlinear_term = (y(1).^2) .* y(2);
    
    derivative_u = a - y(1) + nonlinear_term;
    derivative_v = b - nonlinear_term;

    derivative = [derivative_u, derivative_v]';
end


%This event function is used to terminate the simulation, if it exceeds a
%suitably large threshold.
function [position,isterminal,direction] = terminalEventFcn(threshold, t,y, Z)
    position = max(abs(y(1)))<threshold; %The value that we want to be zero
    isterminal = 1;  % Halt integration 
    direction = 0; 
end
