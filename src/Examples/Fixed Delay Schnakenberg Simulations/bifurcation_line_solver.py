# @brief    This files is designed to be an example of how to produce
#           bifurcation plots for the fixed delay Schnakenberg system.

import numpy as NP;
import scipy.linalg as LA;
import matplotlib.pyplot as plt
from pychebfun import chebfun
from scipy.optimize import minimize_scalar
from scipy.optimize import minimize
from scipy.optimize import root_scalar

from scipy import optimize



def computeBValue(a_value, tau):

    
    
    
    
    u_fixed = lambda a,b: b+a;
    v_fixed = lambda a,b: b/((b+a)**2);
    
    
    
    lhs = lambda u, v, beta: 1 - ((2*u*v*(beta**2))/(u**4 + (beta**2) * (u**2 - 2*u*v)**2))**2
    rhs = lambda u,v, beta: ((beta**2)/(u**4)) * (1+(u**2 - 2*u*v)*2*u*v*(beta**2)/(u**4 + (beta**2)*(u**2 - 2*u*v)**2))**2
    
    
    
    
    
    beta_values = NP.linspace(0, 10, num=10);
    b_test_values = NP.linspace(0.1,3,num=5000);
    
    b_values=[];
    
    
    for beta_index in range(0, len(beta_values)):
        for b_index in range(0, len(b_test_values)):
            
            beta_value = beta_values[beta_index];
            b_value = b_test_values[b_index];

            #objective_function = lambda b: -lhs(u_fixed(a_value, b), v_fixed(a_value, b), beta_value) + rhs(u_fixed(a_value, b), v_fixed(a_value, b), beta_value);

            #res = minimize_scalar(objective_function);
        
            residue =  -lhs(u_fixed(a_value, b_value), v_fixed(a_value, b_value), beta_value) + rhs(u_fixed(a_value, b_value), v_fixed(a_value, b_value), beta_value);


            if(abs(residue)<1e-3):            
                b_values = b_values + [b_value];
        
                #print("b values: "+str(residue));
        
    return b_values;




# Generate a list of a and b values that will be used to produce a bifurcation
# plot.    
a_amt = 50;
a_values = NP.linspace(0.001, 3, num=a_amt);

b_amt=50;
b_values = NP.linspace(0.001, 3, num=b_amt);


# The tau values to be used for each of the bifurcation plots.
#tau_values = NP.linspace(0, 10, num=50);
tau_values =[0,1,2];


# The numeric values that will be assigned for each type of fixed point.
source_colour_value = 5;
sink_colour_value = 3;
saddle_colour_value = 10;
unstable_spiral_colour_value = -5;
stable_spiral_colour_value = -3;



#computeEigenvalues(0.5, 0.02, 0);

# Iterate through each of the delay values and produce a bifurcation plot for
# of them.
for tau_index in range(0, len(tau_values)):

    


    current_tau = tau_values[tau_index];

    # The matrix that will be used to store the type of fixed point that our
    # system contains for each combination of a and b values.
    fixed_point_types = NP.zeros((a_amt,b_amt));

    # Iterate through each combination of a and b values in the a_values and
    # b_values lists and classify the fixed point in each system.
    for i in range(0, a_amt):
        a = a_values[i];
   
        b_values = computeBValue(a, current_tau);
    
        if(len(b_values)>0):
            for j in range(0, len(b_values)):
                b_index = (int)(NP.floor(b_values[j]*50/3))-1;
                
                
                fixed_point_types[i][b_index] = 10;
                    
               

    # Plot types of fixed point in a bifurcatiion diagram.
    XX, YY = NP.meshgrid(a_values, b_values);
    fig = plt.figure();
    plt.title(r'Bifurcation Plot for $\tau = '+str(round(current_tau, 2))+'$ Delay Schnakenberg Kinetics');
    p = plt.imshow(fixed_point_types, extent=[0, max(a_values), 3, 0]);
    plt.colorbar(p);
    plt.xlabel('$a$');
    plt.ylabel('$b$');
    plt.gca().invert_yaxis()
    plt.show();
    #fig.savefig('./Fixed_Delay_Bifurcation_Plot_Figures/fig_'+str(tau_index)+'.png', dpi=300);
    #plt.close(fig);
    