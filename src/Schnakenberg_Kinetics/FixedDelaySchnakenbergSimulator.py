##  @brief      This file defines the FixedDelaySchnakenbergSimulator class.


from ddeint import ddeint
import numpy as NP



##  @brief      This class is designed to produce numerical simulations of
#               trajectories in the fixed delay Schnakenberg ODE system:
#               \f$ \frac{du}{dt}(t) = a - u(t) + u^2(t-\tau)v(t-\tau) \f$
#               \f$ \frac{dv}{dt}(t) = b - u^2(t-\tau)v(t-\tau) \f$, for 
#               a given set of parameters.
class FixedDelaySchnakenbergSimulator:
    
    ##  @brief              The constructor of this class.
    #   @param  a           This must be a value of type double, which
    #                       corresponds to the value of the a coefficient
    #                       which will be used in the governing equation:
    #                       \f$ \frac{du}{dt}(t) = a - u(t) + 
    #                       u^2(t-\tau)v(t-tau) \f$.
    #   @param  b           This must be a value of type double, which
    #                       corresponds to the value of the b coefficient
    #                       which will be used in the govenring equation:
    #                       \f$ \frac{dv}{dt}(t) = b - u^2(t-\tau)v(t-\tau)
    #                       \f$.
    #   @param t_values     This must be an array of type double, whose 
    #                       entries are all >=0. This argument corresponds
    #                       to the values in time at which we which we wish
    #                       to determine the value of u(t) and v(t) in our
    #                       trajectory simulations.
    #   @param tau          This must be a value of type double which is >0
    #                       and corresponds to the time delay which will be
    #                       used in the governing equations.
    def __init__(self, a, b, tau, t_values):
        self.__a = a;
        self.__b = b;
        self.__t_values = t_values;
        self.__tau = tau;
        
      

    ##  @brief      This method produces a simulation of a trajectory in
    #               the fixed time delay Schnakenberg ODE system, for a given
    #               set of initial conditions.    
    #   @param u0   This must be a function which takes a value of type double
    #               in the interval \f$ [-\tau, 0] \f$ and returns a value of
    #               type double. This corresponds to the initial data for the
    #               variable u, that we wish to use for the trajectory
    #               simulation.    
    #   @param u0   This must be a function which takes a value of type double
    #               in the interval \f$ [-\tau, 0] \f$ and returns a value of
    #               type double. This corresponds to the initial data for the
    #               variable v, that we wish to use for the trajectory
    #               simulation.    
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
        initial_data = lambda t: (u0(t),v0(t));
        sol = ddeint(self.__evaluateModel, initial_data, self.__t_values);        
       
        u_sol = sol[:,0];
        v_sol = sol[:,1];
    
        return (u_sol, v_sol);


    ##  @brief      This method is used to compute the derivative of u and 
    #               v, using the fixed delay Schnakenberg kinetics for a
    #               given set of previous values of u and v.
    #   @param  y   This must be a function, which takes a value of type 
    #               double, in the interval \f$[t-\tau, t]\f$ and returns
    #               a one dimensional array of length 2, of type
    #               double. This corresponds to the previous values
    #               attained by the functions u and v, prior to time t.
    #   @param  t   This must be a value of type double, which corresponds
    #               to the time at which we wish to compute the derivative
    #               of u and v under the fixed time delay Schnakenberg
    #               kinetics.
    #   @return     This method returns a one dimensional array of length 2,
    #               of type double, which corresponds to the value of the
    #               derivatives of u and v respectively under the dynamics of
    #               fixed time delay Schnakenberg kinetics.
    def __evaluateModel(self, y, t):
        u, v = y(t);
        u_tau, v_tau = y(t-self.__tau);
        return NP.array([self.__a - u + (u_tau**2)*v_tau, 
                         self.__b - (u_tau**2)*v_tau]);
