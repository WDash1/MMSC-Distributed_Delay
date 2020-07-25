% The parameter settings to be used for simulations of the Schnakenberg
% reaction-diffusion system.
a = 0.2;
b = 1.3;

%The diffusion coefficients for simulations of the PDE model.
Du = 1;
Dv = 100;

%The time delay value we wish to use for our model.
tau=0.60;

%The fixed point for our model.
u_fixed = b + a;
v_fixed = b/((b+a).^2);

%The number of points which we wish to divide the space interval into.
n = 50;

%The maximum and minimum spatial values that we wish to use for
%simulations.
x_max = 50;
x_min = 0;

%The discretised space values that will be used in simulations of our PDE.
x_values = (x_min + ((0:n)/n) .*(x_max-x_min))';

%The space step size for our numerical simulations.
h = (x_max - x_min)/n;

%The initial data at t<=0 for our simulation.
u0 = @(x) u_fixed + 0.001*(sin(sqrt(2)*x) + cos(x));
v0 = @(x) v_fixed + 0.001*(sin(sqrt(2)*x) + cos(x));
initial_data = @(t) [u0(x_values).*(u0(t)) ; v0(x_values).*(v0(t))];

%The time derivative function for the delay equation.
dydt = @(t,y,Z) computeFixedDelaySchnakenbergDerivative(a,b,Du,Dv,n,h,y,Z);

%Simulate the model and retrieve the values of the solution on the
%specified time mesh.
t_mesh = linspace(1,100,100);

eventFunction = @(t,y,Z) terminalEventFcn(10e5,t,y,Z);
options = odeset('RelTol',10e-5,'Events',eventFunction);
sol = dde23(dydt,tau,initial_data,t_mesh,options);
solution_mesh = sol.y';

%Extract the u and v values from the simulation.
u_values = solution_mesh(:,1:n+1);
v_values = solution_mesh(:,n+1+(1:(n+1)));

[XX,YY] = meshgrid(x_values,sol.x);

%Plot the u values from the simulation.
figure(1);
h1=surf(XX,YY, u_values);
view([0  90])
axis square
grid off
set(h1, 'EdgeColor','none')
xlabel('x');
ylabel('t');
title('u values');
colorbar

%Plot the v values from the simulation.
figure(2);
h2=surf(XX,YY, v_values);
view([0  90])
axis square
grid off
set(h2, 'EdgeColor','none')
xlabel('x');
ylabel('t');
title('v values');
colorbar


%This event function is used to terminate the simulation, if it exceeds a
%suitably large threshold.
function [position,isterminal,direction] = terminalEventFcn(threshold, t,y, Z)
    position = max(abs(y(1)))<threshold; %The value that we want to be zero
    isterminal = 1;  % Halt integration 
    direction = 0; 
end