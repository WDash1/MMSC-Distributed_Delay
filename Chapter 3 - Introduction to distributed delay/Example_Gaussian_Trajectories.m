t_mesh = linspace(0,0.5,500);


%The mean and standard deviation of the distribution.
mu = 0.1;
sigma = mu*0.1;

%The number of integrand evaluation points.
N=13;

%The limits of integration.
a = 0;
b = mu*2;


%Compute the Gauss Hermite quadrature weights for our given normal
%distribution, which lie inside the integration limits [a,b].
[evaluation_points, integrand_weights] = computeWeights(a,b,mu,sigma,N)

%Use the resulting weights to approximate the integral of a normal
%distribution over the interval [a,b].
approximation = sum(integrand_weights)
%normcdf = @(x) (1/2) .* (1-erf(-x./(sqrt(2))));
actual_value = (normcdf((b-mu)/(sigma)) - normcdf((a-mu)/(sigma)))



%The derivative function for our example distributed delay equation.
dydt = @(t,y,Z) exampleDerivative(integrand_weights(:)', a,b,mu,sigma,t,y,Z);
    
%Generate an event function to stop the simulation, if it exceeds the
%specified threshold.
eventFunction = @(t,y,Z) terminalEventFcn(10e7,t,y,Z);
options = odeset('RelTol',10e-10,'AbsTol',10e-10,'Events',eventFunction);


%The initial data.
u0 = @(t) exp(t);
v0 = @(t) exp(t);
y0 = @(t) [u0(t),v0(t)];



%If zero is an evaluation point in our delay_times vector, then exclude
%it from the dde23 function call and re-insert it in the underlying
%distributedDelaySchnakenberg function.
if(evaluation_points(1) == 0)  
    solution = dde23(dydt, double(evaluation_points(2:end)), y0, [0,max(t_mesh)], options);
else
    solution = dde23(dydt, double(evaluation_points), y0, [0,max(t_mesh)], options);
end


%Plot the expected solution as well as the resulting approximations.
y_values = deval(solution, t_mesh);
figure
hold on
plot(t_mesh, y_values(1,:), 'DisplayName', 'x(t)');
plot(t_mesh, y_values(2,:), 'DisplayName', 'y(t)');
plot(t_mesh, exp(t_mesh), 'DisplayName', 'exp(t)');
legend
xlabel("t");
hold off

%Plot the error associated with the resulting approximations.
figure
hold on

x_err = abs(y_values(1,:) - exp(t_mesh))
y_err = abs(y_values(2,:) - exp(t_mesh))

plot(t_mesh, abs(y_values(1,:) - exp(t_mesh)));
plot(t_mesh, abs(y_values(2,:) - exp(t_mesh)));
hold off
xlabel("t");
ylabel("Approximation Error");
legend("x error", "y error");


%This function is used to compute the Gauss Hermite weights and evaluation
%points and then filter our the values which lie outside of the domain 
%[a,b]
function [evaluation_points, integrand_weights] = computeWeights(a,b,mu,sigma,N)
    [points, weights] = computeGaussHermiteWeights(mu,sigma,N);
    
    
    valid_indicies = (points > a) .* (points < b);
    
    evaluation_points = points .* valid_indicies;
    integrand_weights = weights .* valid_indicies;
    
    evaluation_points = nonzeros(evaluation_points);
    integrand_weights = nonzeros(integrand_weights);
    
end


%This function is used to compute the derivative of the example distributed
%delay equation.
function derivative = exampleDerivative(weights, a, b, mu, sigma, t,y,Z)
    function_values = Z;
    if(size(Z,2)+1 == size(weights,2))
        function_values = [y, Z];
    end
    
    
    integrand = function_values(2,:).^2;
    integral_approximation = dot(weights, integrand);
    
    
    
    mu1 = (2*t*sigma - mu) / (2*sigma.^2 - 1);
    sigma1 = sigma/(sqrt(1-2*sigma^2));
    
    
    derivative_x = y(2);
    derivative_y = y(1) + integral_approximation - ...
                    (exp(t.^2 - mu.^2/(2.*sigma.^2))./(sigma.*sqrt(2*pi))) .* ...
                    exp(-((2*sigma.^2 - 1)/(2.*sigma)).*(((mu - 2*t*sigma.^2)./(2*sigma.^2 - 1)).^2)) .* ...
                    sqrt(2*pi)*sigma1*(normcdf(b, mu1, sigma1)- normcdf(a,mu1,sigma1));
                    %sqrt(2*pi)*sigma1*(normcdf(b, mu1, sigma1)- normcdf(a,mu1,sigma1));
    
    
    
                
     %derivative_y = y(1) + integral_approximation - y(2).^2 .* ...
     %               exp((((2*(sigma^2) - mu).^2)-mu.^2)/(sigma.^2)) .* ...
     %               (normcdf(b+2*sigma.^2, mu, sigma) - normcdf(a+2*sigma.^2, mu, sigma));
    
     derivative = [double(derivative_x), double(derivative_y)]';
end


function [position,isterminal,direction] = terminalEventFcn(threshold, t,y, Z)
    position = max(abs(y(1)))<threshold; %The value that we want to be zero
    isterminal = 1;  % Halt integration 
    direction = 0; 
end