%The time mesh on which we evaluate trajectories.
t_mesh = linspace(0,40,50000);

%Parameter values
a = 0.05;
b = 0.5;
tau_amt = 4;
tau_values = linspace(0.01,0.04,tau_amt);


%The fixed point for the parameter values.
u_fixed = b + a;
v_fixed = b/((b+a).^2);

%The initial data.
u0 = @(t) u_fixed + 0.1*(sin(sqrt(2)*t) + cos(t));
v0 = @(t) v_fixed + 0.1*(sin(sqrt(2)*t) + cos(t));



[no_delay_u_values, no_delay_v_values] = computeNoDelaySchnakenbergTrajectory(a,b,t_mesh,u0(0),v0(0),10e-12,10e9);
%plot(u_values, v_values, 'DisplayName',"No delay");

delay_u_values = zeros(size(tau_values,2), size(t_mesh,2));
delay_v_values = zeros(size(tau_values,2), size(t_mesh,2));



for tau_index = 1:tau_amt
    tau_value = tau_values(tau_index);
    
    %Plot the fixed delay trajectory.
    sol = computeFixedDelayRBSchnakenbergTrajectory(a,b,tau_value,max(t_mesh), u0,v0,10e-12,10e9);
    current_y_values = deval(sol, t_mesh);
    delay_u_values(tau_index,:) = current_y_values(1,:);
    delay_v_values(tau_index,:) = current_y_values(2,:);
end




%Plot the u-v trajectories
figure('Renderer', 'painters', 'Position', [10 10 320 320], 'Visible', 'on')
set(gca,'FontSize',10)
box on
hold on;
xlabel('{\it u}');
ylabel('{\it v}');
plot(no_delay_u_values, no_delay_v_values, 'DisplayName', 'No Delay');
for tau_index = 1:tau_amt
    plot(delay_u_values(tau_index,:)', delay_v_values(tau_index,:)', 'DisplayName',"\tau="+string(tau_values(tau_index)));
end
plot(u_fixed, v_fixed, '*r', 'DisplayName', 'Fixed Point');
%legend
hold off;
filename = "Output_Images/Fixed_Delay_RB_Model_uv_Trajectory.eps";
print('-depsc', '-tiff', '-r300', '-painters', filename);



%Plot the u trajectories
figure('Renderer', 'painters', 'Position', [10 10 320 320], 'Visible', 'on')
set(gca,'FontSize',10)
box on
hold on;
xlabel('{\it t}');
ylabel('{\it u}');
plot(t_mesh, no_delay_u_values,'DisplayName', 'No Delay');
for tau_index = 1:tau_amt
    plot(t_mesh, delay_u_values(tau_index,:), 'DisplayName',"\tau="+string(tau_values(tau_index)));
end
hold off;
legend
filename = "Output_Images/Fixed_Delay_RB_Model_u_values.eps";
print('-depsc', '-tiff', '-r300', '-painters', filename);


%Plot the v trajectories
figure('Renderer', 'painters', 'Position', [10 10 320 320], 'Visible', 'on')
set(gca,'FontSize',10)
box on
hold on;
xlabel('{\it t}');
ylabel('{\it v}');
plot(t_mesh, no_delay_v_values,'DisplayName', 'No Delay');
for tau_index = 1:tau_amt
    plot(t_mesh, delay_v_values(tau_index,:), 'DisplayName',"\tau="+string(tau_values(tau_index)));
end
hold off;
legend
filename = "Output_Images/Fixed_Delay_RB_Model_v_values.eps";
print('-depsc', '-tiff', '-r300', '-painters', filename);


