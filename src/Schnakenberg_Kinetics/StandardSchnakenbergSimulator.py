##  @brief      This file defines the StandardSchnakenbergSimulator class.


from scipy.integrate import odeint



##  @brief      This class is designed to produce numerical simulations of
#               trajectories in the standard Schnakenberg ODE system:
#               \f$ \frac{du}{dt}(t) = a - u(t) + u^2(t)v(t) \f$
#               \f$ \frac{dv}{dt}(t) = b - u^2(t)v(t) \f$, for a given set of
#               parameters.
class StandardSchnakenbergSimulator:
    
    ##  @brief              The constructor of this class.
    #   @param  a           This must be a value of type double, which
    #                       corresponds to the value of the a coefficient
    #                       which will be used in the governing equation:
    #                       \f$ \frac{du}{dt}(t) = a - u(t) + u^2(t)v(t) \f$.
    #   @param  b           This must be a value of type double, which
    #                       corresponds to the value of the b coefficient
    #                       which will be used in the govenring equation:
    #                       \f$ \frac{dv}{dt}(t) = b - u^2(t)v(t) \f$.
    #   @param t_values     This must be an array of type double, whose 
    #                       entries are all >=0. This argument corresponds
    #                       to the values in time at which we which we wish
    #                       to determine the value of u(t) and v(t) in our
    #                       trajectory simulations.
    #   @param  atol        This is an optional parameter which corresponds
    #                       to the absolute tolerance used by the ODE solver
    #                       for simulating trajectories. For more details,
    #                       please consult the funciton odeint.
    #   @param  rtol        This is an optional parameter which corresponds to
    #                       the absolute tolerance used by the ODE solver for
    #                       simulating trajectories. For more details,
    #                       please consult the function odeint.
    def __init__(self, a, b, t_values, atol=None, rtol=None):
        self.__a = a;
        self.__b = b;
        self.__t_values = t_values;
        self.__abserr = atol;
        self.__relerr = rtol;
        
        self.__f = lambda y, t: (self.__a - y[0] + (y[0]**2) * y[1], 
                                 self.__b - (y[0]**2) * y[1]);

    ##  @brief      This function produces a simulation of a trajectory in
    #               the standard Schnakenberg ODE system, for a given set of 
    #               initial conditions.    
    #   @param u0   This must be a value of type double, which corresponds to
    #               the initial condition for the u variable at time 0, that
    #               we wish to use for the trajectory simulation.    
    #   @param v0   This must be a value of type double, which corresponds to 
    #               the initial condition for the v variable at time 0, that
    #               we wish to use for the trajectory simulation.    
    #   @return     This function returns two arrays of size n, where n
    #               is used to denote the length of the t_values array given
    #               to the constructor of this class. The first of these
    #               arrays corresponds to the values of u attained on the
    #               requested trajectory at the time values specified in the
    #               t_values array given to the constructor of this class.
    #               The second of these arrays corresponds to the value of v
    #               attained on the requested trajectory at the time values
    #               specified in the t_values array given to the constructor
    #               of this class.
    def generateTrajectory(self, u0, v0):
        sol = odeint(self.__f, [u0,v0], self.__t_values,
                     atol=self.__abserr, rtol=self.__relerr);
   
        u_sol = sol[:,0];
        v_sol = sol[:,1];
    
        return (u_sol, v_sol);