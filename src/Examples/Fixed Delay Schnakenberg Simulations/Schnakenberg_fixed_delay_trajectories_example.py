##  @brief      This file acts as an example of how one might produce various
#               plots of trajectories in the fixed delay Schnakenberg 
#               system, by using the FixedDelaySchnakenbergSimulator
#               class.

import numpy as NP;
from pylab import figure, plot, xlabel, ylabel, legend, title, savefig
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

import sys
sys.path.append('../Schnakenberg_Kinetics/')
from FixedDelaySchnakenbergSimulator import FixedDelaySchnakenbergSimulator;


# Parameter values for the kinetics.
a=9
b=2
tau = 1;


# Compute the corresponding fixed point values for the system.
u_fixed = b+a;
v_fixed = b/((b+a)**2);


# The time values at which we wish to compute the values of the trajectories.
t_values = NP.linspace(0, 2, num=300);

# Create a trajectory simulator object with a fixed delay.
trajectory_simulator = FixedDelaySchnakenbergSimulator(a, b, tau, t_values);


# Simulate trajectories for the standard Schnakenberg system, using
# several different initial conditons.
u0_values = [lambda t: t, lambda t: 5*t, lambda t: 0.8*t];
v0_values = [lambda t: 2*t, lambda t: t,lambda t: 8*t];
u_sol_1, v_sol_1 = trajectory_simulator.generateTrajectory(u0_values[0], 
                                                         v0_values[0]);
u_sol_2, v_sol_2 = trajectory_simulator.generateTrajectory(u0_values[1],
                                                         v0_values[1]);
u_sol_3, v_sol_3 = trajectory_simulator.generateTrajectory(u0_values[2],
                                                         v0_values[2]);


# The style, colour and name settings for each line in the subsequent plots.
line_colours = ['r', 'g', 'b', 'y'];
line_styles = ['-', '-', '-', '--'];
line_settings = list(map(lambda x,y: x+y, line_colours, line_styles));
legend_strings = ['$u_0 = t, v_0 = 2t$',
        '$u_0 = 5t, v_0 = t$',
        '$u_0 = 0.8t, v_0 = 8t$',
        '$u_0 = u_*, v_0 = v_*$'];

# Plot the u values over time.
figure(1);
xlabel('$t$');
ylabel('$u$');
plot(t_values, u_sol_1, line_settings[0]);
plot(t_values, u_sol_2, line_settings[1]);
plot(t_values, u_sol_3, line_settings[2]);
plot(t_values, NP.ones((len(t_values),1)) * u_fixed, line_settings[3]);
legend(legend_strings);
title(r'Fixed $\tau=$'+str(tau)+" Delay Schnakenberg $u(t)$ Simuations");


# Plot the v values over time.
figure(2);
xlabel('$t$');
ylabel('$v$');
plot(t_values, v_sol_1, line_settings[0]);
plot(t_values, v_sol_2, line_settings[1]);
plot(t_values, v_sol_3, line_settings[2]);
plot(t_values, NP.ones((len(t_values),1)) * v_fixed, line_settings[3]);
legend(legend_strings);
title(r'Fixed $\tau=$'+str(tau)+' Delay Schnakenberg $v(t)$ Simuations');



# Plot the trajectory curves in 2D.
figure(3);
xlabel('$u$');
ylabel('$v$');
plot(u_sol_1, v_sol_1, line_settings[0]);
plot(u_sol_2, v_sol_2, line_settings[1]);
plot(u_sol_3, v_sol_3, line_settings[2]);
plot([u_fixed], [v_fixed], line_colours[3]+'o');
legend(legend_strings);
title(r'Fixed $\tau=$'+str(tau)+' Delay Schnakenberg Trajectories');



# Plot the trajecotry curves in 3D.
fig = plt.figure(4);
ax = plt.axes(projection='3d');
ax.set_xlabel('$u$');
ax.set_ylabel('$v$');
ax.set_zlabel('$t$');
ax.plot(u_sol_1, v_sol_1, t_values, line_settings[0]);
ax.plot(u_sol_2, v_sol_2, t_values, line_settings[1]);
ax.plot(u_sol_3, v_sol_3, t_values, line_settings[2]);
ax.plot(NP.ones((len(t_values), 1)) * u_fixed,
        NP.ones((len(t_values), 1)) * v_fixed,
        t_values, line_settings[3])
legend(legend_strings);
fig.suptitle(r'Fixed $\tau=$'+str(tau)+' Delay Schnakenberg Trajectories');