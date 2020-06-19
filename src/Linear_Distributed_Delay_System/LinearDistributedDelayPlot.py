import numpy as NP
from pylab import figure, plot, xlabel, ylabel, legend, title, savefig
from DistributedDelaySimulator import DistributedDelaySimulator;


# The number of points we wish to use the Trapezium rule to discretise the
# integral.
n=100;

# The model parameters we wish to use for simulations.
alpha=-2
beta_amt=5;
beta_values = NP.linspace(-11, 0, num=beta_amt);

# The time values at which we wish to compute the values of the trajectories.
t_values = NP.linspace(0, 20, num=300);

# Initial data for the simulations.
y0_values = lambda t: NP.sin(NP.sqrt(2)*t)+NP.cos(t);

# Produce trajectory simulations for each value of beta.
distributed_delay_simulator = DistributedDelaySimulator(n, t_values);
y_sols = list(map((lambda beta_value: distributed_delay_simulator.generateTrajectory(y0_values, alpha, beta_value)), beta_values));

# The style, colour and name settings for each line in the subsequent plots.
line_colours = ['r', 'g', 'b', 'y', 'm'];
line_styles = ['-', '-', '-', '-', '-'];
line_settings = list(map(lambda x,y: x+y, line_colours, line_styles));
legend_strings=[];

# Plot the y values over time, for each set of model parameters.
figure(1);
xlabel('$t$');
ylabel('$y$');
for i in range(0, beta_amt):
    legend_strings.append(r'$ \beta = '+str(beta_values[i])+'$');
    plot(t_values, y_sols[i], line_settings[i]);

legend(legend_strings);
title(r'Linear Distributed delay Trajectory Simuations for $\alpha = '+str(alpha)+'$');

