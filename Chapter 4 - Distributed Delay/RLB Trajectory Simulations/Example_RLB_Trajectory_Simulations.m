%The time mesh on which we evaluate trajectories.
t_mesh = linspace(00,40,50000);

%Parameter values
a = 0.05;
b = 0.5;
%a=1
%b=1
tau = 0.04;

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
[no_delay_u_values, no_delay_v_values] = computeNoDelaySchnakenbergTrajectory(a,b,t_mesh,u0(0),v0(0),10e-18,10e9);


%Plot the fixed delay trajectory.
sol = computeFixedDelayRLBSchnakenbergTrajectory(a,b,tau,max(t_mesh), u0,v0,10e-18,10e9);
current_y_values = deval(sol,t_mesh);
fixed_delay_u_values = current_y_values(1,:);
fixed_delay_v_values = current_y_values(2,:);


sigma_amt = 3;
distributed_delay_u_values = zeros(sigma_amt, size(t_mesh,2));
distributed_delay_v_values = zeros(sigma_amt, size(t_mesh,2));

N=11;
sigma_values = [tau*0.01, tau*0.1, tau*0.2];
%sigma_values = linspace(tau*0.1,tau*0.1,sigma_amt);
%sigma_values = linspace(tau*0.01,tau*0.1,sigma_amt);
for sigma_index = 1:sigma_amt
    current_sigma = sigma_values(sigma_index);


    [delay_values, weights] = computeGaussHermiteWeights(0,2*tau,tau,current_sigma,N)

    sol = computeDistributedDelayRLBSchnakenbergTrajectory(a,b,double(delay_values), double(weights), max(t_mesh), u0,v0,10e-18,10e9);
    current_y_values = deval(sol,t_mesh);
    distributed_delay_u_values(sigma_index, :) = current_y_values(1,:);
    distributed_delay_v_values(sigma_index, :) = current_y_values(2,:);
end



%Plot the u-v trajectories
figure('Renderer', 'painters', 'Position', [10 10 320 320], 'Visible', 'on')
set(gca,'FontSize',10)
box on
hold on;
xlabel('{\it u}');
ylabel('{\it v}');
plot(no_delay_u_values, no_delay_v_values, 'DisplayName', 'No Delay');
plot(fixed_delay_u_values, fixed_delay_v_values, 'DisplayName', '\sigma=0');
for sigma_index = 1:sigma_amt
    plot(distributed_delay_u_values(sigma_index,:)', distributed_delay_v_values(sigma_index,:)', 'DisplayName',"\sigma="+string(sigma_values(sigma_index)));
end
plot(u_fixed, v_fixed, '*r', 'DisplayName', 'Fixed Point');
%legend
hold off;
filename = "Output_Images/Distributed_Delay_RLB_Model_uv_Trajectory.eps";
print('-depsc', '-tiff', '-r300', '-painters', filename);


%Plot the u trajectories
figure('Renderer', 'painters', 'Position', [10 10 320 320], 'Visible', 'on')
set(gca,'FontSize',10)
box on
hold on;
xlabel('{\it t}');
ylabel('{\it u}');
plot(t_mesh, no_delay_u_values, 'DisplayName', 'No Delay');
plot(t_mesh, fixed_delay_u_values, 'DisplayName', '\sigma=0');
for sigma_index = 1:sigma_amt
    plot(t_mesh, distributed_delay_u_values(sigma_index,:)', 'DisplayName',"\sigma="+string(sigma_values(sigma_index)));
end
%legend
hold off;
filename = "Output_Images/Distributed_Delay_LI_Model_u_values.eps";
print('-depsc', '-tiff', '-r300', '-painters', filename);


%Plot the v trajectories
figure('Renderer', 'painters', 'Position', [10 10 320 320], 'Visible', 'on')
set(gca,'FontSize',10)
box on
hold on;
xlabel('{\it t}');
ylabel('{\it v}');
plot(t_mesh, no_delay_v_values, 'DisplayName', 'No Delay');
plot(t_mesh, fixed_delay_v_values, 'DisplayName', '\sigma=0');
for sigma_index = 1:sigma_amt
    plot(t_mesh, distributed_delay_v_values(sigma_index,:)', 'DisplayName',"\sigma="+string(sigma_values(sigma_index)));
end
%legend
hold off;
filename = "Output_Images/Distributed_Delay_RLB_Model_v_values.eps";
print('-depsc', '-tiff', '-r300', '-painters', filename);



