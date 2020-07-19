%The time values for which we wish to simulate the trajectories of the
%system.
t_values = linspace(10,30,500);

%The number of integral discretisation points.
N=10;


%The parameter values we wish to use.
epsilon = 0;
delta_amt = 50;
gamma_amt = 50;
delta_values = linspace(-5*pi^2,20*pi^2,delta_amt);
gamma_values = linspace(-50*pi^2,20*pi^2,gamma_amt);

%The initial condtions we wish to use.
y0 = @(t) [sin(sqrt(2)*t)+cos(t); sin(sqrt(2)*t)+cos(t)];

%Run simulations for each combination of parameters.
Ms = zeros(gamma_amt, delta_amt);
parfor i=1:delta_amt
    current_delta = delta_values(i);
    for j=1:gamma_amt
        current_gamma = gamma_values(j);
        [solx,soly] = compute_trajectory_simulation(epsilon, current_gamma, current_delta, N, t_values, y0);
        if(max(abs(soly(1,:)))>1)
            Ms(j,i) = 6;
        else
            Ms(j,i) = 0;
        end
    end
end


%Put the matrix into a figure.
figure('Renderer', 'painters', 'Position', [10 10 300 300], 'Visible', 'on')
colormap(jet);
imagesc(linspace(min(delta_values./(pi^2)),max(delta_values./(pi^2)),delta_amt), linspace(min(gamma_values./(pi^2)),max(gamma_values./(pi^2)),gamma_amt), Ms);
caxis([0,10]);
%colorbar;
    
set(gca,'YDir','normal')
xlabel('$\frac{\delta}{\pi^2}$', 'interpreter','latex');
ylabel('$\frac{\gamma}{\pi^2}$', 'interpreter','latex');
set(gca,'FontSize',10)
    
%x_axis_labels = min(alpha_values):1:max(alpha_values);
%y_axis_labels = min(beta_values):1:max(beta_values);
    
set(gca,'XTick',min(delta_values./(pi^2)):5:max(delta_values./(pi^2)));
set(gca,'xticklabel',num2str(get(gca,'xtick')','%.0f'))
    
set(gca,'YTick',min(gamma_values./(pi^2)):10:max(gamma_values./(pi^2)));
set(gca,'yticklabel',num2str(get(gca,'ytick')','%.0f'))
    
%    print('-r300','-dpng','Fig'+string(tau_index)+'.png');
print('Linear_Example', '-dpng', '-r300');
    




%This function produces a trajectory simulation, for a given set of
%parameters.
function [x,y] = compute_trajectory_simulation(epsilon, gamma, delta, N, t_values, y0)
    K = @(s) pi/2*sin(pi*s);

    delay_times = linspace(0,1,N);
    dydt = @(t,y,Z) model(K, delay_times, epsilon, gamma, delta, N, t,y,Z);
    delay_times2 = delay_times(2:end);

    sol = dde23(dydt, delay_times2, y0, t_values);
    
    x = sol.x;
    y = sol.y;
end

%This function returns the derivatie of the model.
function derivative = model(kernel_function, delay_times, epsilon, gamma, delta, N, t,y,Z)
    y_delay_values1 = [y(1), Z(1,:)];
    y_delay_values2 = [y(2), Z(2,:)];

  
    kernel_values = kernel_function(delay_times);
    
    sum_values = sum(kernel_values .* y_delay_values1);
    derivative = [y(2); -(delta + epsilon * cos(4*pi*t)) * y(1) + gamma/N*(-(y(1)+y_delay_values1(end))/2 + sum_values)];
end