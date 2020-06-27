git diff --name-only --diff-filter=U
import numpy as NP
from pylab import figure, plot, xlabel, ylabel, legend, title, savefig
import matplotlib.pyplot as plt
from DistributedDelaySimulator import DistributedDelaySimulator;


# The number of points we wish to use the Trapezium rule to discretise the
# integral.
n=20;

# The model parameters we wish to use for simulations.
alpha_amt = 50
alpha_values = NP.linspace(-10, 2, num=alpha_amt);

beta_amt = 50;
beta_values = NP.linspace(-20, 10, num=beta_amt);

# The time values at which we wish to compute the values of the trajectories.
t_values = NP.linspace(10, 20, num=300);

# Initial data for the simulations.
y0_values = lambda t: NP.sin(NP.sqrt(2)*t)+NP.cos(t);

# Produce trajectory simulations for each value of beta.
distributed_delay_simulator = DistributedDelaySimulator(n, t_values);




stability_matrix = NP.zeros((beta_amt, alpha_amt));

for i in range(0, alpha_amt):
    for j in range(0, beta_amt):
        alpha_value = alpha_values[i];
        beta_value = beta_values[j];
        
        y_sol = distributed_delay_simulator.generateTrajectory(y0_values, 
                                                               alpha_value, 
                                                               beta_value);
        
        
        if(max(abs(y_sol))>1.0):
            stability_matrix[j][i] = 6;
        else:
            stability_matrix[j][i] = 0;


# Plot types of fixed point in a bifurcatiion diagram.
XX, YY = NP.meshgrid(alpha_values, beta_values);
fig,ax = plt.subplots(1,1)

plt.title('Bifurcation Plot for Linear Distributed Delay System');
p = plt.imshow(stability_matrix, 
              extent=[min(alpha_values), max(alpha_values), max(beta_values), min(beta_values)], 
              aspect = (max(alpha_values)-min(alpha_values))/(max(beta_values)-min(beta_values)),
              cmap=plt.cm.get_cmap('jet'));

plt.clim(0,10)
fig, ax = plt.subplots(1,1)
plt.title('Bifurcation Plot for Linear Distributed Delay System');
p = ax.imshow(stability_matrix, 
              extent=[min(alpha_values), max(alpha_values), max(beta_values), min(beta_values)], 
              aspect = (max(alpha_values)-min(alpha_values))/(max(beta_values)-min(beta_values)));

#plt.colorbar(p);
plt.xlabel(r'$\alpha$');
plt.ylabel(r'$\beta$');
plt.gca().invert_yaxis()
fig.savefig('fig2.png', dpi=300);
plt.show();

