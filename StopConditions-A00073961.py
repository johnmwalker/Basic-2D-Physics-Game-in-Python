# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 12:03:20 2016

@author: walk0106
"""
import Simulation
import math

def lessThanX(x):
    def lessThan(y):
        return y < x
    return lessThan
        
def greaterThanX(x):
    def greaterThan(y):
        return y > x
    return greaterThan
        
def numOrbitsSetup(stop,desiredOrbits):
    '''Sets up a stop condition that tells a simulation to stop after a desired
    number of orbits.
    
    Parameters
    -----
    stop : int
        The number of steps at which to trigger a stop in order to prevent an
        infinite loop.
    desiredOrbits: int
        The number of orbits through which you want to simulate.
    
    Returns
    -----
    numOrbits : function
        The stop condition function.
    '''
    def numOrbits(sim):
        '''A stop condition that takes in an instance of a simulation class and
        tells it to stop after a desired number of orbits.
        
        Parameters
        -----
        sim : instance of a simulation class
            The simulation class whose stopping point is being found.
            
        Returns
        -----
        False until the condition is triggered and then it returns True
        '''
        if len(sim.t)>stop:
            sim.numberOrbits = desiredOrbits
            print "Didn't reach desired number of orbits!"
            return True
            
        else: return sim.numberOrbits >= desiredOrbits
    return numOrbits
    
def afterTimeSetup(time=[None],stop=False,periodMultiplier=None,a=None):
    '''Sets up a stop condition that tells a simulation to stop after a 
    specific amount of time or after a desired ratio of the orbiting object's 
    period.
    
    Parameters
    -----
    time : float
        The time in seconds at which the simulation should stop
    periodMultiplier : float
        The multiplier for the period; 0.5 will cause the stop condition to 
        trigger at half of a period (apoapsis/periapsis).
    a : float
        The semimajoraxis of the orbiting body. This is used to calculate the 
        period.
    stop : int
        The number of steps at which to trigger a stop in order to prevent an
        infinite loop. This funcitonality is disabled by default.
    
    Returns
    -----
    afterTime : function
        The stop condition function.
    '''

            
    def afterTime(sim):
        '''A stop condition that takes in an instance of a simulation class and
        tells it to stop after a ratio of the period of the object's orbit.
        
        Parameters
        -----
        sim : instance of a simulation class
            The simulation class whose stopping point is being found.
            
        Returns
        -----
        False until the condition is triggered and then it returns True
        '''
        if periodMultiplier!=None:
            p=math.sqrt(a**3*4*math.pi**2/(sim.M*6.7408*10**-11))
            time[0]=periodMultiplier*p
        if stop != False:
            if len(sim.t)>stop:
                print "Didn't reach desired number of orbits!"
                return True

        else:  
            
            if sim.solver.stepsize > time-sim.t[-1]:
                sim.solver.stepsize = time-sim.t[-1]
                
            return time <= sim.t[-1]
            
    return afterTime

def numStepsSetup(stop,steps):
    '''Sets up a stop condition that tells a simulation to stop after a desired
    number of steps. Note that the condition will trigger when the number of 
    steps the simulation has completed is equal to the number of steps entered.
    
    Parameters
    -----
    steps: int
        The number of steps at which you want to stop the simulation.
    stop : int
        The number of steps at which to trigger a stop in order to prevent an
        infinite loop.
    
    Returns
    -----
    numSteps : function
        The stop condition function.
    '''
    def numSteps(sim):
        '''A stop condition that takes in an instance of a simulation class and
        tells it to stop after a desired number of steps.
        
        Parameters
        -----
        sim : instance of a simulation class
            The simulation class whose stopping point is being found.
            
        Returns
        -----
        False until the condition is triggered and then it returns True
        '''
        if len(sim.t)>stop:
            print "Didn't reach desired number of orbits!"
            return True
            
        else: return len(sim.t)==steps
        
    return numSteps

    
def lessThanZero(x):
    '''Tell your simulation method when to stop! This checks when input < 0.
    
    Parameters
    -----
    x : float
        Input variable to be checked.
    
    Returns
    -----
    False until condition is met, which then returns True
    '''
    return x < 0
    
def numberOfOrbits(bodies,desiredOrbits,numOrbits):
    '''Tell your simulation method when to stop! This takes in orbital info
    and a desired number of orbits and will return True after the simulation
    has completed the desired number of orbits.
    
    Parameters
    -----
    bodies : array of GravBody objects
        All of the bodies at each timestep in a simulation.
    desiredOrbits : int
        The number of orbits to simulate through in total.
    numOrbits : int
        The number of orbits that the simulation has gone through thus far.
    '''

    if len(bodies)<3:
        return numOrbits
    elif len(bodies)>50000:
        print "Didn't reach desired number of orbits!"
        return desiredOrbits
    else:
        originalPos = bodies[0].position
        stopSpan=abs(bodies[1].position-originalPos)
        currentPos = bodies[-1].position
        diff = abs(currentPos-originalPos)
        
#        print 'diff ' + repr(diff.r)
#        print 'stopSpan ' + repr(stopSpan.r)

        if diff.r <= stopSpan.r/2.:
            numOrbits=numOrbits+1
#            print numOrbits
        return numOrbits

