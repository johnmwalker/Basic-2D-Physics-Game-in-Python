# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 14:51:17 2016

@author: John
"""

class Euler(object):
    '''A solver that uses Euler's method.
    
    Attributes
    -----
    stepsize : float
        The size between steps for Euler's method.
    diffEq : class method
        A method containing information about a differential equation.
    '''
    def __init__(self, stepsize,diffEq):
        self.stepsize = stepsize
        self.diffEq = diffEq
        
    def advance(self,F,x):
        '''Advances the simulation based on Euler's method one step forward.
        
        Parameters
        -----
        F : float
            The current y value used in the differential equation.
        x : float
            The current x value used in the differential equation.
            
        Returns
        -----
        FNext : float
            The next y value after using Euler's method.
        xNext : float
            The next x value after one step.
        '''
        FNext = F + self.diffEq(F,x)*self.stepsize
        return FNext, x+self.stepsize
        
class RK2(object):
    '''A solver that uses second order Runge-Kutta (RK2).
    
    Attributes
    -----
    stepsize : float
        The size between steps for RK2.
    diffEq : class method
        A method containing information about a differential equation.
    '''
    def __init__(self, stepsize,diffEq):
        self.stepsize = stepsize
        self.diffEq = diffEq
        
    def advance(self,F,x):
        '''Advances the simulation based on a second order Runge-Kutta (RK2)
        one step forward.
        
        Parameters
        -----
        F : float
            The current y value used in the differential equation.
        x : float
            The current x value used in the differential equation.
            
        Returns
        -----
        FNext : float
            The next y value after using RK2.
        xNext : float
            The next x value after one step.
        '''
        step =  self.stepsize
        k1 = self.diffEq(F,x)
        k2 = step*self.diffEq(F + k1*step *0.5 , x + 0.5*step)
        FNext = F + self.diffEq(F + self.diffEq(F,x)*self.stepsize *0.5,x)*self.stepsize 
        return FNext, x+self.stepsize
        
class RK4(object):
    '''A solver that uses fourth order Runge-Kutta (RK4).
    
    Attributes
    -----
    stepsize : float
        The size between steps for RK4.
    diffEq : class method
        A method containing information about a differential equation.
    '''
    def __init__(self, stepsize,diffEq):
        self.stepsize = stepsize
        self.diffEq = diffEq
        
    def advance(self,F,x):
        '''Advances the simulation based on a fourth order Runge-Kutta (RK4)
        one step forward.
        
        Parameters
        -----
        F : float
            The current y value used in the differential equation.
        x : float
            The current x value used in the differential equation.
            
        Returns
        -----
        FNext : float
            The next y value after using RK4.
        xNext : float
            The next x value after one step.
        '''
        step =  self.stepsize
        k1 = step*self.diffEq(F,x)
        k2 = step*self.diffEq(F + k1*0.5 , x + 0.5*step)
        k3 = step*self.diffEq(F + k2*0.5 , x + 0.5*step)
        k4 = step*self.diffEq(F + k3 , x + step)
        FNext = F + (1/6.)*k1 + (1/3.)*k2 + (1/3.)*k3 + (1/6.)*k4
        return FNext, x+self.stepsize

