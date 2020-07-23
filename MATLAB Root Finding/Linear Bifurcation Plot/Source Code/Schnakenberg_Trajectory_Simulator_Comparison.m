%The time mesh on which we evaluate trajectories.
t_mesh = linspace(0,20,500);

%Parameter values
a = 0.01;
b = 0.5;
tau = 0.06;

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
[u_values, v_values] = computeNoDelaySchnakenbergTrajectory(a,b,t_values,u0(0),v0(0),10e-9,10e9);
plot(u_values, v_values, 'DisplayName',"No delay");

%Plot the fixed delay trajectory.
sol = computeFixedDelaySchnakenbergTrajectory(a,b,tau,max(t_values), u0,v0,10e-9,10e9);
current_y_values = deval(sol, t_mesh);
u_values = current_y_values(1,:);
v_values = current_y_values(2,:);
plot(u_values, v_values, 'DisplayName',"Fixed delay");
    



%The number of points to use in the quadrature rule.
N=39;
[legendre_zeros, legendre_weights] = computeGaussLegendreWeights(0, 2*tau, N);

%The sigma values we wish to simulate.
sigma_amt = 4;
sigma_values = linspace(0.001,0.04,sigma_amt);

%Plot the distributed delay trajectories.
for index=1:sigma_amt
    sigma = sigma_values(index);
    sol = computeDistributedDelaySchnakenbergTrajectory(legendre_zeros, ...
                                                        legendre_weights, ...
                                                        a,b,sigma,...
                                                        max(t_values),...
                                                        u0,v0,...
                                                        10e-9,10e9);
    y_values = deval(sol, t_mesh);
    u_values = y_values(1,:);
    v_values = y_values(2,:);
    
    
    plot(u_values,v_values, 'DisplayName',"\sigma = "+string(round(sigma,10)));
    
end
legend
%print('Trajcetory Comparison Plot', '-dpng', '-r300');


