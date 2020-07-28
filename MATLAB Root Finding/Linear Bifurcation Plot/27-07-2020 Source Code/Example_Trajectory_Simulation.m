%The time mesh on which we evaluate trajectories.
t_mesh = linspace(0,95,5000);

%Parameter values
a = 0.2;
b = 1.3;
tau = 0.75;

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
    




