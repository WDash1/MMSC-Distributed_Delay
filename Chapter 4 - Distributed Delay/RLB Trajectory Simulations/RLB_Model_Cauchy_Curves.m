t_values = linspace(0,0.2,500);

u0 = @(t) sin(sqrt(2)*t)+cos(t);
v0 = @(t) sin(sqrt(2)*t)+cos(t);

a = 3;
b = 3;

mu = 0.04;
sigma = 0.008;

n_max = 10;
n_values = 1:n_max;

trapezium_rule_cauchy_sequence = zeros(size(n_values));
simpsons_rule_cauchy_sequence = zeros(size(n_values));
simpsons_38_rule_cauchy_sequence = zeros(size(n_values));
gauss_legendre_cauchy_sequence = zeros(size(n_values));
gauss_hermite_cauchy_sequence = zeros(size(n_values));

trapezium_current_trajectory = zeros(2,size(t_values,2))';
simpsons_current_trajectory = zeros(2,size(t_values,2))';
simpsons_38_current_trajectory = zeros(2,size(t_values,2))';
gauss_legendre_current_trajectory = zeros(2,size(t_values,2))';
gauss_hermite_current_trajectory = zeros(2,size(t_values,2))';



    
for n = n_values
    spmd(4)
        if labindex == 1
            [trapezium_rule_zeros, trapezium_rule_weights] = computeTrapeziumRuleWeights(0, 2*mu, 6*n+1);
            [~, trapezium_soly] = compute_trajectory_simulation(trapezium_rule_zeros, trapezium_rule_weights,a, b, mu, sigma, t_values, u0, v0);
            if(n>1)
                trapezium_rule_cauchy_sequence(n-1) = max(norm(trapezium_soly - trapezium_current_trajectory));
            end
            trapezium_current_trajectory = trapezium_soly;
        elseif labindex == 2
            [simpsons_rule_zeros, simpsons_rule_weights] = computeCompositeSimpsonRuleWeights(0, 2*mu, 6*n+1);
            [~, simpsons_soly] = compute_trajectory_simulation(simpsons_rule_zeros, simpsons_rule_weights,a, b, mu, sigma, t_values, u0, v0);
            if(n>1)
                simpsons_rule_cauchy_sequence(n-1) = max(norm(simpsons_soly - simpsons_current_trajectory));
            end
            simpsons_current_trajectory = simpsons_soly;
        elseif labindex == 5
            [simpsons_38_rule_zeros, simpsons_38_rule_weights] = computeCompositeSimpson38RuleWeights(0, 2*mu, 6*n+1);
            [~, simpsons_38_soly] = compute_trajectory_simulation(simpsons_38_rule_zeros, simpsons_38_rule_weights,a, b, mu, sigma, t_values, u0, v0);
            if(n>1)
                simpsons_38_rule_cauchy_sequence(n-1) = max(norm(simpsons_38_soly - simpsons_38_current_trajectory));
            end
            simpsons_38_current_trajectory = simpsons_38_soly;
        elseif labindex == 4
            [gauss_legendre_zeros, gauss_legendre_weights] = computeGaussLegendreWeights(0, 2*mu, 6*n+1);
            [~, legendre_soly] = compute_trajectory_simulation(gauss_legendre_zeros, gauss_legendre_weights,a, b, mu, sigma, t_values, u0, v0);
            if(n>1)
                gauss_legendre_cauchy_sequence(n-1) = max(norm(legendre_soly - gauss_legendre_current_trajectory));
            end
            gauss_legendre_current_trajectory = legendre_soly;
        elseif labindex == 3
            [gauss_hermite_zeros, gauss_hermite_weights] = computeGaussLegendreWeights(0, 2*mu, 6*n+1);
            sol = computeDistributedDelayLISchnakenbergTrajectory(...
                                                a, ...
                                                b, ...
                                                gauss_hermite_zeros, ...
                                                gauss_hermite_weights, ...
                                                max(t_values), ...
                                                u0, ...
                                                v0, ...
                                                10e-10, ...
                                                10e10);
            hermite_soly=deval(sol, t_values);
            if(n>1)
                gauss_hermite_cauchy_sequence(n-1) = max(norm(hermite_soly - gauss_hermite_current_trajectory));
            end
            gauss_hermite_current_trajectory = hermite_soly;
        end
        
    end
end


figure('Renderer', 'painters', 'Position', [10 10 500 500], 'Visible', 'on')
hold on;
box on;
xlim([min(6.*n_values+1), max(6.*(n_max-1)+1)]);
loglog(6.*n_values+1, trapezium_rule_cauchy_sequence{1}, '-o', 'DisplayName', 'Trapezium Rule');
loglog(6.*n_values+1, simpsons_rule_cauchy_sequence{2}, '-o', 'DisplayName', "Simpson's Rule");
%loglog(6.*n_values+1, simpsons_38_rule_cauchy_sequence{5}, '-o', 'DisplayName', "Simpson's 3/8 Rule");
loglog(6.*n_values+1, gauss_legendre_cauchy_sequence{4}, '-o', 'DisplayName', "Gauss Legendre");
loglog(6.*n_values+1, gauss_hermite_cauchy_sequence{3}, '-o', 'DisplayName', "Gauss Hermite");
xlabel('{\it N}', 'FontSize', 20);
legend
set(gca,'XScale', 'log', 'YScale', 'log')
ax = gca;
ax.FontSize = 20; 
filename = "LI_Distributed_Delay_Cauchy_Curve_a="+string(a)+"_b="+string(b)+"_tau="+string(mu)+"_sigma="+string(sigma);
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