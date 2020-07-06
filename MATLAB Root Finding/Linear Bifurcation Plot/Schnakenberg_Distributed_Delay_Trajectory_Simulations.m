%The time values for which we wish to simulate the trajectories of the
%system.
t_values = linspace(0,30,100);

%The number of integral discretisation points.
N=61;

%The parameter values that we wish to use for our simulation.
a_amt = 3;
b_amt = 3;
a_values = linspace(0.01,2,a_amt);
b_values = linspace(0.01,10,b_amt);

%The initial data.
y0 = @(t) [sin(sqrt(2)*t)+cos(t), sin(sqrt(2)*t)+cos(t)];

%Run trajectory simulations over the parameter space.
Ms = zeros(b_amt, a_amt);
parfor i=1:a_amt
    current_a = a_values(i);
    for j=1:b_amt
        current_b = b_values(j);
        [solx,soly] = compute_trajectory_simulation(current_a, current_b,1, N, t_values, y0);
        if(max(abs(soly))>1)
            Ms(j,i) = 6;
        else
            Ms(j,i) = 0;
        end
    end
end

%Plot the stability matrix as a figure.
figure('Renderer', 'painters', 'Position', [10 10 300 300], 'Visible', 'on')
colormap(jet);
imagesc(linspace(min(alpha_values),max(alpha_values),alpha_amt), linspace(min(beta_values),max(beta_values),beta_amt), Ms);
caxis([0,10]);
    
set(gca,'YDir','normal')
xlabel('{ \alpha}');
ylabel('{ \beta}');
set(gca,'FontSize',10)
    
%x_axis_labels = min(alpha_values):1:max(alpha_values);
%y_axis_labels = min(beta_values):1:max(beta_values);
    
set(gca,'XTick',0:0.5:max(a_values));
set(gca,'xticklabel',num2str(get(gca,'xtick')','%.1f'))
    
set(gca,'YTick',0:2:max(b_values));
set(gca,'yticklabel',num2str(get(gca,'ytick')','%.0f'))
    
%    print('-r300','-dpng','Fig'+string(tau_index)+'.png');
print('Linear_Example', '-dpng', '-r300');
    



function [position,isterminal,direction] = terminalEventFcn(t,y, Z)
    position = max(abs(y(1)))<10^3; % The value that we want to be zero
    isterminal = 1;  % Halt integration 
    direction = 0; 
end




%This function produces a trajectory simulation, for a given set of
%parameters.
function [x,y] = compute_trajectory_simulation(a, b, K, N, t_values, y0)
    delay_times = linspace(0,1,N);
    delay_times = delay_times(2:end);
    dydt = @(t,y,Z) simpsons_38_rule_model(a, b, K, N, t,y,Z);
    
    options = odeset('RelTol',1e-7,'Events',@terminalEventFcn);
    sol = dde23(dydt, delay_times, y0, t_values, options);
    
    x = sol.x;
    y = sol.y;
end

%This function uses the trapezium rule to approximate the derivative in the
%system, for a given set of parameters.
function derivative = trapezium_rule_model(a, b, K, N, t,y,Z)
    
    integrand = (Z(1,:).^2).*Z(2,:);

    
    sum_value = ((y(1).^2*y(2)) + (y(end).^2*y(end)))/2 + (sum(integrand)-integrand(end));
    
    derivative_u = a - y(1) + 1/N*sum_value;
    derivative_v = b - 1/N*sum_value;
    derivative = [derivative_u, derivative_v]';
end

%This function uses the composite simpsons rule to approximate the
%derivative in the system, for a given set of parameters.
function derivative = simpsons_rule_model(a, b, K,N, t,y,Z)
    even_indicies = linspace(3, N-2, ((N-1)/2)-1);
    odd_indicies = linspace(2, N-1, ((N-1)/2));
    
    
    integrand = (Z(1,:).^2).*Z(2,:);

    sum_value = (y(1).^2*y(2))+2*sum(integrand(even_indicies)) + 4*sum(integrand(odd_indicies))+integrand(end);
    derivative_u = a - y(1) + 1/(3*N)*sum_value;
    derivative_v = b - 1/(3*N)*sum_value;
    derivative = [derivative_u, derivative_v]';
end

%This function uses the composite simpsons 38 rule to approximate the
%derivative in the system, for a given set of parameters.
function derivative = simpsons_38_rule_model(a, b, K, N, t,y,Z)
    third_indicies = linspace(4,N-3, ((N-1)/3)-1);
    
    first_third_indicies = linspace(3,N-1, ((N-1)/3));
    second_third_indicies = linspace(3,N-1, ((N-1)/3));
    
    integrand = (Z(1,:).^2).*Z(2,:);
    
    
    sum_value = (y(1).^2*y(2)) + 3*(sum(integrand(first_third_indicies))+sum(integrand(second_third_indicies)))...
                +2*sum(integrand(third_indicies)) + integrand(end);
            
            
            
    derivative_u = a - y(1) + 3/(8*N)*sum_value;
    derivative_v = b - 3/(8*N)*sum_value;
    derivative = [derivative_u, derivative_v]';
end