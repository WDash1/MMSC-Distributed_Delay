##  @brief      This file defines the FixedDelaySchnakenbergSimulator class.


from ddeint import ddeint
import numpy as NP



##  @brief      This class is designed to produce numerical simulations of
#               trajectories in the linear distributed delay system:
#               \f$ \frac{yu}{dt}(t) = \alpha y(t) + \beta + \int_{0}^{1} y(t-s) ds\f$,
#               for a given set of parameters. This is implement by
#               using the trapezium rule to approximate the integral
#               in this expression to reduce it to a fixed delay
#               differential equation.
class DistributedDelaySimulator:
    
    ##  @brief              The constructor of this class.
    #   @param  n           This must be an integer which is >0 and dictates
    #                       the number of uniformly spaced points in the
    #                       inverval [0,1] at which the integrand is evaluated
    #                       in trajectory simulations, when using the 
    #                       trapezium rule.
    #   @param t_values     This must be an array of type double, whose 
    #                       entries are all >=0. This argument corresponds
    #                       to the values in time at which we which we wish
    #                       to determine the value of y(t) in our
    #                       trajectory simulations.
    def __init__(self, n, t_values):
        self.__n = n;
        self.__t_values = t_values;
        self.__evaluation_points = NP.linspace(0, -1, num = self.__n);

      
    ##  @brief          This method produces a simulation of a trajectory in
    #                   the distributed delay linear system:
    #                   \f$ \frac{yu}{dt}(t) = \alpha y(t) + \beta + \int_{0}^{1} y(t-s) ds\f$,
    #                   for a given set of parameters and initial conditions.
    #   @param alpha    This must be a value of type double, which corresponds
    #                   to the value of the coefficient alpha in the governing
    #                   equation.
    #   @param beta    This must be a value of type double, which corresponds
    #                   to the value of the coefficient beta in the governing
    #                   equation.
    #   @param y0       This must be a function which takes a value of type
    #                   double in the interval \f$ [-1, 0] \f$ and returns a
    #                   value of type double. This corresponds to the initial
    #                   data for the varaible y, that we wish to use for the
    #                   trajectory simulation.
    #   @return         This function returns an arrays of size m, where m
    #                   is used to denote the length of the t_values array
    #                   given to the constructor of this class. The value in
    #                   this array correspond to the values of y attained on
    #                   the requested trajectory at the time values specified
    #                   in the t_values array given to the constructor of this
    #                   class.
    def generateTrajectory(self, y0, alpha, beta):
        initial_data = lambda t: y0(t);
        sol = ddeint(self.__evaluateModel, initial_data, self.__t_values,
                     fargs=(alpha, beta));        
    
        return sol;


    ##  @brief          This method is used to approximate the derivative of y
    #                   in the numerically discretised version of the
    #                   governing equation, using the trapezium rule.
    #   @param  y       This must be a function, which takes a value of type 
    #                   double, in the interval \f$[t-1, t]\f$ and returns
    #                   a value of type double.
    #                   This corresponds to the previous values
    #                   attained by the function y, prior to time t.
    #   @param  t       This must be a value of type double, which corresponds
    #                   to the time at which we wish to compute the derivative
    #                   of the dependent variable y.
    #   @param alpha    This must be a value of type double, which corresponds
    #                   to the value of the coefficient alpha in the governing
    #                   equation.
    #   @param beta    This must be a value of type double, which corresponds
    #                   to the value of the coefficient beta in the governing
    #                   equation.
    #   @return         This method returns a value of type double, which 
    #                   approximates the value of the derivative of y, with 
    #                   respect to time, at time t.
    def __evaluateModel(self, y, t, alpha, beta):
        time_values = self.__evaluation_points + t;
        y_values = list(map(y, time_values));
    

        derivative = alpha*y(t) + beta/self.__n*(-(y(t) + y(t-1))/2 + sum(y_values));

        return derivative;



