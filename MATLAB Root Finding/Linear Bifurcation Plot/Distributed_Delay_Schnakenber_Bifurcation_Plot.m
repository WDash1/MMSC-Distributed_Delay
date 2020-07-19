%The time values for which we wish to simulate the trajectories of the
%system.
t_values = linspace(10,20,10);

tau = 0.03;
sigma = 0.01;


a_amt = 2;
a_values = [0.1,2];
b_amt = 2;
b_values = [0.5,2];

N=30;
[legendre_zeros, legendre_weights] = computeGaussLegendreWeights(0, 2*tau, N);



Ms = zeros(b_amt, a_amt);
parfor i=1:a_amt
    current_a = a_values(i);
    for j=1:b_amt
        current_b = b_values(j);
        
        u_fixed = current_b + current_a;
        v_fixed = current_b/((current_b+current_a).^2);
    
        u0 = @(t) u_fixed + 0.001*(sin(sqrt(2)*t) + cos(t));
        v0 = @(t) v_fixed + 0.001*(sin(sqrt(2)*t) + cos(t));
        y0 = @(t) [u0(t), v0(t)];

        [u_sol, v_sol] = compute_trajectory_simulation(current_a,current_b,sigma,double(legendre_zeros), double(legendre_weights),t_values,y0);
     
        u_deviation = max(abs(u_sol-u_fixed));
        v_deviation = max(abs(v_sol-v_fixed));
        total_deviation = max(u_deviation, v_deviation);

        if(total_deviation<0.01)
            Ms(j,i) = 0;
        else
            Ms(j,i) = 6;
        end
    end
end



%Plot the stability matrix as a figure.
figure('Renderer', 'painters', 'Position', [10 10 300 300], 'Visible', 'on')
colormap(jet);
imagesc(linspace(0,max(a_values),a_amt), linspace(0,max(b_values),b_amt), Ms);
caxis([0,10]);
    
set(gca,'YDir','normal')
xlabel('{\it a}');
ylabel('{\it b}');
set(gca,'FontSize',10)
    
%x_axis_labels = min(alpha_values):1:max(alpha_values);
%y_axis_labels = min(beta_values):1:max(beta_values);
    
set(gca,'XTick',0:0.5:max(a_values));
set(gca,'xticklabel',num2str(get(gca,'xtick')','%.1f'))
    
set(gca,'YTick',0:0.5:max(b_values));
set(gca,'yticklabel',num2str(get(gca,'ytick')','%.1f'))
    
%    print('-r300','-dpng','Fig'+string(tau_index)+'.png');
print('Schnakenberg_Distributed_Delay_Bifurcation_Plot', '-dpng', '-r300');
    



%This event function is used to terminate the simulation, if it exceeds a
%suitably large threshold.
function [position,isterminal,direction] = terminalEventFcn(t,y, Z)
    position = max(abs(y(1)))<10^9; % The value that we want to be zero
    isterminal = 1;  % Halt integration 
    direction = 0; 
end


%This function produces a trajectory simulation, for a given set of
%parameters.
function [u_solution, v_solution] = compute_trajectory_simulation(a,b,sigma,delay_times, weights, t_values, y0)

    t_max = max(t_values);
    tau = max(delay_times)/2;

    phi = @(x) exp(-(x.^2) ./ 2).* 1/(sqrt(2*pi));
    Phi = @(x) (1/2).*(erf(x./sqrt(2)));

    
    K = @(s) 1/sigma .* phi((s - tau) ./ sigma) ./ (Phi(tau / sigma) - Phi(-tau/sigma));
    
    kernel_values = K(delay_times);
    
    dydt = @(t,y,Z) custom_rule_model(delay_times(:)', weights(:)', a,b,kernel_values(:)',t,y,Z);
    
    options = odeset('RelTol',1e-10,'Events',@terminalEventFcn);
    
    if(delay_times(1) == 0)  
        sol = dde23(dydt, delay_times(2:end), y0, [0,t_max], options);
    else
        sol = dde23(dydt, delay_times, y0, [0,t_max], options);
    end
    
    y_solution = deval(sol, t_values);
    u_solution = y_solution(1,:);
    v_solution = y_solution(2,:);
end


%This function uses Gauss-Legendre quadrature to approximate the
%derivative in the system, for a given set of parameters.
function derivative = custom_rule_model(delay_values, weights, a, b, kernel_values,t, y, Z)


    function_values = Z;
    if(size(Z,2)+1 == size(delay_values,2))
        function_values = [y, Z];
    end
    
    integrand = ((function_values(1,:)).^2).* function_values(2,:) .* kernel_values;
    integral_approximation = dot(weights, integrand);
    
    derivative_u = a - y(1) + integral_approximation;
    derivative_v = b - integral_approximation;
    derivative = [derivative_u, derivative_v]';
end



