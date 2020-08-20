%The time values for which we wish to simulate the trajectories of the
%system.
t_values = linspace(10,25,500);

%The number of integral discretisation points.
N=11;

%The parameter values that we wish to use for our simulation.
alpha_amt = 500;
beta_amt = 500;
%alpha_values = linspace(-2.5,-1.75,alpha_amt);
%beta_values = linspace(3.25,3.75,beta_amt);
alpha_values = linspace(-10,2,alpha_amt);
beta_values = linspace(-20,10,beta_amt);

%The initial data.
y0 = @(t) sin(sqrt(2)*t)+cos(t);


[gauss_legendre_zeros, gauss_legendre_weights] = computeGaussLegendreWeights(0, 1, N);



%Run trajectory simulations over the parameter space.
Ms = zeros(beta_amt, alpha_amt);
parfor i=1:alpha_amt
    current_alpha = alpha_values(i);
    for j=1:beta_amt
        current_beta = beta_values(j);
        [solx,soly] = compute_trajectory_simulation(gauss_legendre_zeros, gauss_legendre_weights, current_alpha, current_beta, t_values, y0);
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
    
set(gca,'XTick',min(alpha_values):2:max(alpha_values));
set(gca,'xticklabel',num2str(get(gca,'xtick')','%.0f'))
    
set(gca,'YTick',min(beta_values):5:max(beta_values));
set(gca,'yticklabel',num2str(get(gca,'ytick')','%.0f'))
    
%    print('-r300','-dpng','Fig'+string(tau_index)+'.png');
print('Linear_Example', '-dpng', '-r300');
    

function [x,y] = compute_trajectory_simulation(delay_times,weights, alpha, beta, t_values, y0)
    delays = delay_times;
    if(delay_times(1) == 0)
        delays = delay_times(2:end);
    end
    dydt = @(t,y,Z) linearDerivativeExample(weights,alpha, beta, t,y,Z);
    options = odeset('RelTol',10e-10,'AbsTol',10e-10, 'Events', @terminalEventFcn);
    sol = dde23(dydt, delays, y0, t_values,options);
    
    x = sol.x;
    if(size(sol.ye,1) == 0)
        y = deval(sol, t_values);
    else
        y = sol.y;
    end
    
end


function [position,isterminal,direction] = terminalEventFcn(t,y, Z)
    position = max(abs(y(1)))<10^3; % The value that we want to be zero
    isterminal = 1;  % Halt integration 
    direction = 0; 
end





function derivative = linearDerivativeExample(weights, alpha, beta, t,y,Z)
    function_values = Z;
    if(size(Z,2)+1 == size(weights,2))
        function_values = [y, Z];
    end
    integrand = function_values;
    integral_approximation = dot(weights, integrand);
    derivative = alpha*y + beta*integral_approximation;
end
