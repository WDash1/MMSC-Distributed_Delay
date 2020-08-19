

a_amt = 200;
b_amt = 200;

a_values = linspace(0.01, 5, a_amt);
b_values = linspace(0.01, 5, b_amt);

tau_amt  = 30;
tau_values = linspace(0, 0.3, tau_amt);

%tau_values = [0.1];    
Ms = zeros(a_amt, b_amt);
for tau_index = 1:tau_amt
    current_tau = tau_values(tau_index);
    parfor i=1:a_amt
        current_a = a_values(i);
        for j=1:b_amt
            current_b = b_values(j);
            Ms(j,i) = 6.0*(DispersRel(current_tau, current_a, current_b)>0);
        end
    end



    %Need to sort out axis labels/plotting here! Blue is no instability, yellow
    %is unstable.
    %figure('Visible', 'on')
    figure('Renderer', 'painters', 'Position', [10 10 300 300], 'Visible', 'on')
    
   
    colormap(jet);
    
    %image('XData',[0 max(As)],'YData',[0 max(Bs)],'CData',Ms) ;
    imagesc(linspace(0,max(a_values),a_amt), linspace(0,max(b_values),b_amt), Ms);
    caxis([0,10]);
    %colorbar;
    
    set(gca,'YDir','normal')
    xlabel('{\it a}');
    ylabel('{\it b}');
    set(gca,'FontSize',15)
    
    x_axis_labels = 0:1:max(a_values);
    y_axis_labels = 0:1:max(b_values);
    
    set(gca,'XTick',0:1:max(a_values));
    set(gca,'xticklabel',num2str(get(gca,'xtick')','%.1f'))
    
    set(gca,'YTick',0:1:max(b_values));
    set(gca,'yticklabel',num2str(get(gca,'ytick')','%.1f'))
    
    title(strcat('\tau = ', num2str(current_tau, '%.2f')));
%    print('-r300','-dpng','Fig'+string(tau_index)+'.png');
    print('Output_Images/Fig'+string(tau_index), '-dpng', '-r300');
end
