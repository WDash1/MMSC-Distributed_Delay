import numpy as NP;

from ddeint import ddeint

from pylab import figure, plot, xlabel, ylabel, legend, title, savefig
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d








def fullModel(y, t, tau, a, b):
    
    u, v = y(t);
    u_tau, v_tau = y(t-tau);
    return NP.array([a - u + (u_tau**2)*v_tau, 
                     b - (u_tau**2)*v_tau]);

def linearisedModel(y, t, tau, a, b):
    
    u_fixed = b+a;
    v_fixed = b/((b+a)**2);
    
    u, v = y(t);
    u_tau, v_tau = y(t-tau);
    return NP.array([- u + 2*u_fixed*v_fixed*u_tau + (u_fixed**2)*v_tau, 
                     - 2*u_fixed*v_fixed*u_tau - (u_fixed**2)*v_tau]);





def simulateTrajectory(a, b, tau, t_values, u0, v0):
    initial_data = lambda t: (u0(t),v0(t));
    sol = ddeint(linearisedModel, initial_data, t_values, fargs=(tau, a, b));        
       
    u_sol = sol[:,0];
    v_sol = sol[:,1];
    
    return (u_sol, v_sol);

    




def determineStability(a,b,tau):
    t_values=NP.linspace(10,20,num=5);
    trajectory_simulator = FixedDelaySchnakenbergSimulator(a, b, tau, t_values);


    u_fixed = b+a;
    v_fixed = b/((b+a)**2);


    u0_values = lambda t: NP.sin(NP.sqrt(2)*t)+NP.cos(t);
    v0_values = lambda t: NP.sin(NP.sqrt(2)*t)+NP.cos(t);

    u_sol_1, v_sol_1 = simulateTrajectory(a, b, tau, t_values, u0_values, 
                                                         v0_values);
    
     
    
    
    max_u_deviation = max(abs(u_sol_1));
    max_v_deviation = max(abs(v_sol_1));
    
    
    if(max_u_deviation>1 or max_v_deviation>1):
        return False;
    
    return True;







# Generate a list of a and b values that will be used to produce a bifurcation
# plot.    
a_amt = 100;
a_values = NP.linspace(0.001, 3, num=a_amt);
b_amt=100;
b_values = NP.linspace(0.001, 3, num=b_amt);



# The matrix that will be used to store the type of fixed point that our
# system contains for each combination of a and b values.
fixed_point_types = NP.zeros((a_amt,b_amt));


tau_values = NP.linspace(0, 10, num=50);


# The numeric values that will be assigned for each type of fixed point.
source_colour_value = 3;
sink_colour_value = 5;
saddle_colour_value = 10;
unstable_spiral_colour_value = -5;
stable_spiral_colour_value = -3;


for tau_index in range(1, len(tau_values)):

    current_tau = tau_values[tau_index];

    # The matrix that will be used to store the type of fixed point that our
    # system contains for each combination of a and b values.
    fixed_point_types = NP.zeros((a_amt,b_amt));

    # Iterate through each combination of a and b values in the a_values and
    # b_values lists and classify the fixed point in each system.
    for i in range(0, a_amt):
        for j in range(0, b_amt):
            a = a_values[i];
            b = b_values[j];

            
            if(determineStability(a,b,current_tau)):
                fixed_point_types[j][i] = sink_colour_value;
            else:
                fixed_point_types[j][i] = source_colour_value;
                    
               

    # Plot types of fixed point in a bifurcatiion diagram.
    XX, YY = NP.meshgrid(a_values, b_values);
    fig = plt.figure();
    plt.title(r'Bifurcation Plot for $\tau = '+str(round(current_tau, 2))+'$ Delay Schnakenberg Kinetics');
    p = plt.imshow(fixed_point_types, extent=[0, max(a_values), max(b_values), 0]);
    #plt.colorbar(p);
    plt.xlabel('$a$');
    plt.ylabel('$b$');
    plt.gca().invert_yaxis()
    plt.show();
    fig.savefig('./Fixed_Delay_Bifurcation_Plot_Figures2/fig_'+str(tau_index)+'.png', dpi=300);
    plt.close(fig);
        
    
    
    
    