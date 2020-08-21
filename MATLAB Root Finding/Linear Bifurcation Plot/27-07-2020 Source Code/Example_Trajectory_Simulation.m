%The time mesh on which we evaluate trajectories.
t_mesh = linspace(00,100,50000);

%Parameter values
a = 0.3;
b = 0.5;
%a=1
%b=1
tau = 0.89;

%The fixed point for the parameter values.
u_fixed = b + a;
v_fixed = b/((b+a).^2);

%The initial data.
u0 = @(t) u_fixed + 0.1.*(sin(sqrt(2)*t) + cos(t));
v0 = @(t) v_fixed + 0.1.*(sin(sqrt(2)*t) + cos(t));


%Plot the no delay trajectory.
figure
hold on;
xlabel('{\it t}');
ylabel('{\it u}');
plot(v_fixed,u_fixed,'r*', 'DisplayName','Fixed Point')
[u_values1, v_values1] = computeNoDelaySchnakenbergTrajectory(a,b,t_mesh,u0(0),v0(0),10e-18,10e9);
%plot(u_values1,v_values1, 'DisplayName',"No delay");
plot(t_mesh,u_values1, 'DisplayName',"No delay");


%Plot the fixed delay trajectory.
sol = computeFixedDelaySchnakenbergTrajectory(a,b,tau,max(t_mesh), u0,v0,10e-18,10e9);
%current_y_values = deval(sol, t_mesh);
current_y_values = deval(sol,t_mesh)
u_values = current_y_values(1,:);
v_values = current_y_values(2,:);
%plot(u_values,v_values, 'DisplayName',"Fixed delay");
plot(t_mesh,u_values, 'DisplayName',"Fixed delay");


sigma_amt = 5;
N=21;
sigma_values = linspace(tau*0.02,tau*0.2,sigma_amt);
for sigma_index = 1:sigma_amt
    current_sigma = sigma_values(sigma_index);


    [delay_values, weights] = computeGaussHermiteWeights(0,2*tau,tau,current_sigma,N);
    sol = computeDistributedDelaySchnakenbergTrajectory(a,b,double(delay_values), double(weights), max(t_mesh), u0,v0,10e-18,10e9);
    current_y_values = deval(sol,t_mesh)
    u_values = current_y_values(1,:);
    v_values = current_y_values(2,:);
    %plot(u_values,v_values, 'DisplayName',"\sigma = "+string(current_sigma));
    plot(t_mesh,u_values, 'DisplayName',"\sigma = "+string(current_sigma));
end

legend




