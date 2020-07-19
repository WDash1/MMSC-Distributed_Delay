t_values = linspace(10,40,500);

y0 = @(t) sin(sqrt(2)*t)+cos(t);


alpha_values = -2.125;
beta_amt = 10;
beta_values = linspace(3.5625, 3.6875, beta_amt);



for i = 1:beta_amt
    alpha_value = alpha_values;
    beta_value = beta_values(i)
    [solx, soly] = compute_trajectory_simulation(alpha_value, beta_value, N, t_values, y0);
    hold on
    plot(solx, soly);
end
legend({strcat('\beta = ',num2str(beta_values(1))), ...
    strcat('\beta = ',num2str(beta_values(2))), ...
    strcat('\beta = ',num2str(beta_values(3))), ...
    strcat('\beta = ',num2str(beta_values(4))), ...
    strcat('\beta = ',num2str(beta_values(5))), ... 
    strcat('\beta = ',num2str(beta_values(6))), ...
    strcat('\beta = ',num2str(beta_values(7))), ...
    strcat('\beta = ',num2str(beta_values(8))), ...
    strcat('\beta = ',num2str(beta_values(9))), ...
    strcat('\beta = ',num2str(beta_values(10)))}, 'Location', 'southwest');

beta_amt = 10;
N=20;



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