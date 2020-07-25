% The parameter settings to be used for simulations of the Schnakenberg
% reaction-diffusion system.
a = 0.2;
b = 1.3;

%The diffusion coefficients for simulations of the PDE model.
Du = 1;
Dv = 119;

%The fixed point for our model.
u_fixed = b + a;
v_fixed = b/((b+a).^2);

%The number of points which we wish to divide the space interval into.
n = 100;

%The maximum and minimum spatial values that we wish to use for
%simulations.
x_max = 100;
x_min = 0;

%The discretised space values that will be used in simulations of our PDE.
x_values = (x_min + ((0:n)/n) .*(x_max-x_min))';

%The space step size for our numerical simulations.
h = (x_max - x_min)/n;

%The initial data at t=0 for our simulation.
u0 = @(x) u_fixed + 0.001*(sin(sqrt(2)*x) + cos(x));
v0 = @(x) v_fixed + 0.001*(sin(sqrt(2)*x) + cos(x));
initial_data = [u0(x_values) ; v0(x_values)];

%The derivative function for our model.
dydt = @(t,y) computeNoDelaySchnakenbergDerivative(a,b,Du,Dv,n,h,y);

%The time mesh for our model.
t_mesh = linspace(1,100,100);

%Simulate our PDE and extract the u and v values.
[time_mesh,solution_mesh] = ode23(dydt,t_mesh,initial_data);
u_values = solution_mesh(:,1:n+1);
v_values = solution_mesh(:,n+1+(1:(n+1)));
[XX,YY] = meshgrid(x_values,t_mesh);

%Plot the u values from the simulation.
figure(1);
h1=surf(XX,YY, u_values);
view([0  90])
axis square
grid off
set(h1, 'EdgeColor','none')
colorbar

%Plot the v values from the simulation.
figure(2);
h2=surf(XX,YY, v_values);
view([0  90])
axis square
grid off
set(h2, 'EdgeColor','none')
colorbar

