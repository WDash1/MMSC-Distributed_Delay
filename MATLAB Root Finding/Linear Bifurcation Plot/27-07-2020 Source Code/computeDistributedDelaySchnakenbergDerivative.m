%@brief     This function is used to compute the derivative for the
%           reaction diffusion delay PDE model when using fixed delay
%           Shnakenberg reaction kinetics.
%@param a   This must be a scalar which is >0 and dictates the value of the
%           a parameter that we wish to use in the Schnakenberg kinetics.
%@param b   This must be a scalar which is >0 and dictates the value of the
%           b parameter that we wish to use in the Schnakenberg kinetics.
%@param Du  This must be a scalar which is >0 and corresponds to the
%           diffusion coefficient for the u variable in the model.
%@param Dv  This must be a scalar which is >0 and corresponds to the
%           diffusion coefficient for the v variable in the model.
%@param n   This must be an integer which is >0 and dictates the number of
%           points that we wish to segement our 1D space interval into
%           for the discretisation of our model.
%@param h   This must be a scalar which corresponds to the size of each
%           space step in the discretisation of our model.
%
%@param y   This must be a column vector of size 2(n+1), where the first
%           n+1 entries correspond to the current value of u(t) at each of
%           the spatial descretisation points, ordered from the smallest
%           to the largest spatial position. Similarly, the second n+1
%           entries correspond to the current value of v(t) at each of the
%           spatial discretisation points, ordered from the smallest to the
%           largest spatial position.
%@param Z   This must be a column vector of size 2(n+1), which is ordered
%           the same fasion as the y column vector, but corresponds to the
%           current value of the u(t-tau) and v(t-tau) variables.
%@return    This function returns a column vector of size 2(n+1) that
%           dictates the value of the corresponding time derivative for
%           each of the values in the y vector.
function derivative = computeDistributedDelaySchnakenbergDerivative(a,b,Du,Dv, weights,n,h,y,Z)
    function_values = Z;
    if(size(Z,2)+1 == size(weights,2))
        function_values = [y, Z];
    end
    u_delay_values = function_values(1:(n+1), :);
    v_delay_values = function_values(n+1+(1:(n+1)), :);
    integrand = (u_delay_values.^2).*v_delay_values;
    
    %integral_approximation = dot(weights, integrand);
    integral_approximation = integrand*weights;
    

    %The fixed delay Schnakenberg delay differential equation function.
    schnakenberg_u = @(u,v) a - u + -2.*(u.^2.*v) + 3.*integral_approximation;
    schnakenberg_v = @(u,v) b - u.^2.*v;
    
    %The derivative matrix for the spatial dimension, using the central
    %difference approximation.
    second_derivative_matrix = -2.*diag(ones(1,n+1)) + diag(ones(1,n),1) + diag(ones(1,n),-1);
    second_derivative_matrix(1,:) = zeros(1,n+1);
    second_derivative_matrix(n+1,:) = zeros(1,n+1);
    
    %Apply Neumann boundary conditions to the derivative matrix.
    top_row = zeros(1,n+1);
    top_row(1) = -2;
    top_row(2) = 2;

    bottom_row = zeros(1,n+1);
    bottom_row(n) = -2;
    bottom_row(n-1) = 2;

    second_derivative_matrix(1,:) = top_row;
    second_derivative_matrix(n+1,:) = bottom_row;

    %Compute the spatial derivative matrix for the u and v variables.
    u_derivative_matrix = (Du./(h.^2)) .* second_derivative_matrix;
    v_derivative_matrix = (Dv./(h.^2)) .* second_derivative_matrix;

    
    %Add the delay differential equation functions to the spatial
    %derivative matrix, to produce a time derivative function.
    dudt = @(u,v,u_tau,v_tau) u_derivative_matrix * u + schnakenberg_u(u,v);
    dvdt = @(u,v,u_tau,v_tau) v_derivative_matrix * v + schnakenberg_v(u,v);

    u_values = y(1:(n+1));
    v_values =  y(n+1+(1:(n+1)));
   

    derivative = [dudt(u_values, v_values); dvdt(u_values, v_values)];
end