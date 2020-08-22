t_values = linspace(0,10,500);

u0 = @(t) sin(sqrt(2)*t)+cos(t);
v0 = @(t) sin(sqrt(2)*t)+cos(t);

a = 0.2;
b = 1.3;

mu = 5;
%sigma = mu*0.1;

n_amt = 10;
n_values = 5.*(1:n_amt)+1;





sigma_amt = 3;
sigma_values = [mu*0.01, mu*0.1, mu*0.2];

gauss_hermite_cauchy_sequence = zeros(sigma_amt, n_amt);

gauss_hermite_current_trajectory = zeros(2,size(t_values,2), sigma_amt);


for sigma_index = 1:sigma_amt
    
    for n_index = 1:n_amt
        n = n_values(n_index);
        sigma = sigma_values(sigma_index);
        [gauss_hermite_zeros, gauss_hermite_weights] = computeGaussHermiteWeights(0, 2*mu, mu, sigma, n);
        sol = computeDistributedDelayLISchnakenbergTrajectory(...
                                                a, ...
                                                b, ...
                                                double(gauss_hermite_zeros), ...
                                                double(gauss_hermite_weights), ...
                                                max(t_values), ...
                                                u0, ...
                                                v0, ...
                                                10e-10, ...
                                                10e10);
       

        hermite_soly = deval(sol, t_values);
        if(n_index>1)
            gauss_hermite_cauchy_sequence(sigma_index,n_index-1) = max(norm(hermite_soly - gauss_hermite_current_trajectory(:,:,sigma_index)));
        end
        gauss_hermite_current_trajectory(:,:,sigma_index) = hermite_soly;
        
    end
end


figure('Renderer', 'painters', 'Position', [10 10 500 500], 'Visible', 'on')
hold on;
box on;
xlim([min(n_values), max(n_values(1:end-1))]);
for sigma_index = 1:sigma_amt
    loglog(n_values, gauss_hermite_cauchy_sequence(sigma_index,:), '-o', 'DisplayName', '\sigma='+string(sigma_values(sigma_index)));
end
xlabel('{\it N}', 'FontSize', 20);
%legend
set(gca,'XScale', 'log', 'YScale', 'log')
ax = gca;
ax.FontSize = 20; 
filename = "LI_Distributed_Delay_Cauchy_Curve23_a="+string(a)+"_b="+string(b)+"_tau="+string(mu)+"_sigma="+string(sigma);
print('-depsc', '-tiff', '-r300', '-painters', filename+".eps");

function [x,y] = compute_trajectory_simulation(delay_times,weights, a,b,mu,sigma,t_values, u0,v0)
    dist = makedist('Normal', mu,sigma);
    trunc = truncate(dist, 0, 2*mu);
    truncated_gaussian_coefficients = pdf(trunc, delay_times);
    weights2 = weights.*truncated_gaussian_coefficients;

    sol = computeDistributedDelayLISchnakenbergTrajectory(...
                                                a, ...
                                                b, ...
                                                delay_times, ...
                                                weights2, ...
                                                max(t_values), ...
                                                u0, ...
                                                v0, ...
                                                10e-10, ...
                                                10e10);

    x = sol.x;
    y = deval(sol, t_values);
    
end