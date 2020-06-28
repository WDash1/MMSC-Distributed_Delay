t_values = linspace(10,30,500);


alpha_amt = 10;
beta_amt = 10;
N=20;

alpha_values = linspace(-10,2,alpha_amt);
beta_values = linspace(-20,10,beta_amt);

y0 = @(t) sin(sqrt(2)*t)+cos(t);


Ms = zeros(beta_amt, alpha_amt);
parfor i=1:alpha_amt
    current_alpha = alpha_values(i);
    for j=1:beta_amt
        current_beta = beta_values(j);
        [solx,soly] = compute_trajectory_simulation(current_alpha, current_beta, N, t_values, y0);
        if(max(abs(soly))>1)
            Ms(j,i) = 6;
        else
            Ms(j,i) = 0;
        end
    end
end



%figure('Visible', 'on')
figure('Renderer', 'painters', 'Position', [10 10 300 300], 'Visible', 'on')
    
   
colormap(jet);
    
%image('XData',[0 max(As)],'YData',[0 max(Bs)],'CData',Ms) ;
imagesc(linspace(min(alpha_values),max(alpha_values),alpha_amt), linspace(min(beta_values),max(beta_values),beta_amt), Ms);
caxis([0,10]);
%colorbar;
    
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
    





function [x,y] = compute_trajectory_simulation(alpha, beta, N, t_values, y0)
    delay_times = linspace(0,1,N);
    delay_times = delay_times(2:end);
    dydt = @(t,y,Z) model(alpha, beta, N, t,y,Z);
    
    sol = dde23(dydt, delay_times, y0, t_values);
    
    x = sol.x;
    y = sol.y;
end


function derivative = model(alpha, beta, N, t,y,Z)
    sum_values = sum(Z);
    derivative = alpha*y + beta/N*(-(y+Z(end,end))/2 + sum_values);
end