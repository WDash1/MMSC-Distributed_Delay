

a_amt = 10;
b_amt = 10;

a_values = linspace(0.01, 5, a_amt);
b_values = linspace(0.01, 5, b_amt);

tau_amt  = 1;
%tau_values = linspace(0, 1, tau_amt);
tau_values = [0];    
Ms = zeros(M);
for tau_index = 1:tau_amt
    current_tau = tau_values(tau_index);
    parfor i=1:M
        current_a = a_values(i);
        for j=1:M
            current_b = b_values(j);
            Ms(i,j) = 4.0 - 3.25*(DispersRel(current_tau, current_a, current_b)>0);
        end
    end



    %Need to sort out axis labels/plotting here! Blue is no instability, yellow
    %is unstable.
    %figure('Visible', 'on')
    figure('Renderer', 'painters', 'Position', [10 10 800 600], 'Visible', 'on')
    colormap(hsv(128))
    %image('XData',[0 max(As)],'YData',[0 max(Bs)],'CData',Ms) ;
    imagesc(linspace(0,max(a_values),a_amt), linspace(0,max(b_values),b_amt), Ms);
    caxis([0 5]);
    %colorbar;
    
    set(gca,'YDir','normal')
    xlabel('{\it a}');
    ylabel('{\it b}');
    title('Bifurcation Plot for \tau = '+string(current_tau)+' Delay Schnakenberg Kinetics');
%    print('-r300','-dpng','Fig'+string(tau_index)+'.png');
    print('Output_Images/Fig'+string(tau_index), '-dpng', '-r300');
end
