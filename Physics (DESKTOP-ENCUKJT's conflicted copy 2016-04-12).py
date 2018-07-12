# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 14:50:12 2016

@author: John
"""
import Search
import numpy as np
import matplotlib.pyplot as plt
import Body
import Solver
import math
import Simulation
import Vector

class NewtonCooling(object):
    '''Contains methods for use in calculating cooling curves.
    
    Attributes
    -----
    solver : object from Solver class
        Determines the type of approximation you want. Setting to None will
        calculate the exact cooling curve.
    Ta : float
        The ambient temperature for the cooling system.
    k : float
        The cooling constant.
    '''
    def __init__(self,solver,Ta,k):
        self.Ta = Ta
        self.k = k
        if solver != None:
            self.solver = solver
            solver.diffEq = self.diffEq
        
    def diffEq(self, T, t):
        ''' Standard differential equation for cooling, dT/dt=k(Ta-T). Ta comes
        from the class object.
        
        Parameters
        -----
        T : float or int
            The temperature at which you want to find the slope of the cooling curve. (Kelvin)
        k : float or int, default is 0.01
            The cooling constant.

        Returns
        -----
        The result of k(Ta-T).
        '''
        return (self.Ta - T)*self.k
    
    def advance(self,body,time):
        '''Advances the given approximation forward by a step.
        
        Parameters
        -----
        body : a ThermalBody object
            The body whose cooling curve is being calculated.
        time : float
            The time at which the approximation is being advanced.
            
        Returns
        -----
        body : a ThermalBody object
            The object whose temperature has been changed.
        time : float
            The new time after having advanced a step.
        '''
        temperature, time = self.solver.advance(body.temperature,time)
        body.temperature = temperature
        return body, time
        
    def cooling(self,T0, Ta, k, step):
        '''Calculates a cooling curve with a given initial temperature, ambient
        temperature, and cooling constant. Will return an array of temperature 
        values in Kelvin over time from the initial temperature to the ambient 
        temperature.
        
        Parameters
        -----
        T0 : float
            The initial temperature in Kelvin.
        Ta : float
            The ambient temperature in Kelvin.
        k : float
            The cooling constant.
        step : float
            The step (or accuracy) for producing the resulting array.
            
        Returns
        -----
        T : array
            The array of temperature values in Kelvin.
        '''
        e=math.e    
        
        tFinal=-math.log(0.1)/k    # Calculates when the substance is 90% cooled.
        t=np.arange(0,tFinal,step) # Generates the base array on which to do the function.
        T=Ta + (T0-Ta) * e**(-k*t) # Does the final calculation.
        return T
        
class UniformGravity(object):
    '''Contains methods for use in calculating positions and velocities of
    bodies under the influence of gravity.

    Attributes
    -----
    solver : instance of Solver class
        Determines the type of approximation you want.
    a : instance of Vector class
        The acceleration of the system.
    c : float
        The coefficient due to drag. Defaults to zero.
    '''
    def __init__(self,solver,a=Vector.Vector(x=0,y=0,z=-9.81),c=0):
        self.a = a
        self.solver = solver
        self.c = c
        solver.diffEq = self.diffEq
        
    def diffEq(self,f,t):
        ''' Standard differential equations for motion.
        
        Parameters
        -----
        f : numpy array
            A numpy array containing two instances of a Vector class, first for
            velocity, second for position.

        Returns
        -----
        A numpy array containing the velocity and position Vectors as two lists
        within the array.
        '''
        velocity = f[0]
        position = f[1]
        
        vx = velocity[0]
        vy = velocity[1]
        vz = velocity[2]
        
        dvxdt=self.a.x - self.c*vx**2
        dvydt=self.a.y - self.c*vy**2
        dvzdt=self.a.z
    
        dxdt=vx
        dydt=vy
        dzdt=vz
        
        dvdt = [dvxdt,dvydt,dvzdt]
        dpdt = [dxdt,dydt,dzdt]
        
        return np.array([dvdt,dpdt])
    
    def advance(self,body,time):
        '''Advances the given approximation forward by a step.
        
        Parameters
        -----
        body : an instance of GravBody
            The body whose positions and velocities are being calculated.
        time : float
            The time at which the approximation is being advanced.
            
        Returns
        -----
        body : an instance of GravBody
            The object whose positions and velocities are being calculated.
        time : float
            The new time after having advanced a step.
        '''
        x = body.position.x
        y = body.position.y
        z = body.position.z
#        print z
        
        vx=body.velocity.x
        vy=body.velocity.y
        vz=body.velocity.z
        
        pos = np.array([x,y,z])
        vel = np.array([vx,vy,vz])

        f, time = self.solver.advance(np.array([vel,pos]),time)
            
#        print f.shape
#        print f
        
        body.position.x = f[1][0]
        body.position.y = f[1][1]
        body.position.z = f[1][2]
        
        body.velocity.x = f[0][0]
        body.velocity.y = f[0][1]
        body.velocity.z = f[0][2]
        
        return body, time
        
class InverseSquareGravity(object):
    '''Contains methods for use in calculating positions and velocities of
    bodies under the influence of gravity accurately modeled as decreasing as
    an inverse function of r^2.

    Attributes
    -----
    solver : object from Solver class
        Determines the type of approximation you want.
    M : float
        The mass of the base planet or body on which the object is being thrown.
    '''
    def __init__(self,solver,M,a=[0,0,-9.81]):
#        G = -6.67408*10**(-11) # m^3/(kgs^2), the Gravitational Constant
        self.a = a
#        self.M = M
        self.solver = solver
        solver.diffEq = self.diffEq
        
    def diffEq(self,f,t):
        ''' Standard differential equations for motion in a changing
        gravitational field.
        
        Parameters
        -----
        f : numpy array
            A numpy array containing velocity and position values (in that order).

        Returns
        -----
        A numpy array containing the velocity and position values (in that order).
        '''
        velocity = f[0]
        position = f[1]
        
        vx = velocity[0]
        vy = velocity[1]
        vz = velocity[2]
        
        dvxdt=self.a.x
        dvydt=self.a.y
        dvzdt=self.a.z
    
        dxdt=vx
        dydt=vy
        dzdt=vz
        
        dvdt = [dvxdt,dvydt,dvzdt]
        dpdt = [dxdt,dydt,dzdt]
        
        return np.array([dvdt,dpdt])
    
    def advance(self,body,time):
        '''Advances the given approximation forward by a step.
        
        Parameters
        -----
        body : an instance of GravBody
            The body whose positions and velocities are being calculated.
        time : float
            The time at which the approximation is being advanced.
            
        Returns
        -----
        body : an instance of GravBody
            The object whose positions and velocities are being calculated.
        time : float
            The new time after having advanced a step.
        '''
        x = body.position.x
        y = body.position.y
        z = body.position.z
        
        dvxdt=self.a.x
        dvydt=self.a.y
        dvzdt=self.a.z
        
        pos = [x,y,z]
        vel = [dvxdt,dvydt,dvzdt]

        f, time = self.solver.advance(np.array([vel,pos]),time)
            
        body.position.x = f[1][0]
        body.position.y = f[1][1]
        body.position.z = f[1][2]
        
        body.velocity.x = f[0][0]
        body.velocity.y = f[0][1]
        body.velocity.z = f[0][2]
        
        return body, time
        
class CentralGravity(object):
    '''Contains methods for use in calculating positions and velocities of
    bodies under the influence of a central gravity object.

    Attributes
    -----
    solver : object from Solver class
        Determines the type of approximation you want.
    M : float
        The mass of the central gravity object.
    '''
    def __init__(self,solver,M):
        self.G = -6.67408*10**(-11) # m^3/(kgs^2), the Gravitational Constant

        self.M = M
        self.solver = solver
        solver.diffEq = self.diffEq
        
    def diffEq(self,F,t):
        ''' Standard differential equations for motion.
        
        Parameters
        -----
        F : numpy array
            A numpy array containing two instances of a Vector class, first for
            velocity, second for position.

        Returns
        -----
        A numpy array containing the velocity and position Vectors as two lists
        within the array.
        '''
        G = self.G
        M = self.M
        
        wholeLottaStuff=[]
        
        for f in range(len(F)):
#            print 'Body [vel,pos] is ' + repr(F[f])
            velocity = F[f][0]
            pos = F[f][1]
            r = math.sqrt(pos[0]**2+pos[1]**2+pos[2]**2)
            
            vx = velocity[0]
            vy = velocity[1]
            vz = velocity[2]
            
            dvxdt=(G*M/r**3)*pos[0]
            dvydt=(G*M/r**3)*pos[1]
            dvzdt=(G*M/r**3)*pos[2]
        
            dxdt=vx
            dydt=vy
            dzdt=vz
    
            dvdt = [dvxdt,dvydt,dvzdt]
            dpdt = [dxdt,dydt,dzdt]
            
            wholeLottaStuff.append([dvdt,dpdt])
        
        return np.array(wholeLottaStuff)
    
    def advance(self,bodies,time):
        '''Advances the given approximation forward by a step.
        
        Parameters
        -----
        body : an instance of GravBody
            The body whose positions and velocities are being calculated.
        time : float
            The time at which the approximation is being advanced.
            
        Returns
        -----
        body : an instance of GravBody
            The object whose positions and velocities are being calculated.
        time : float
            The new time after having advanced a step.
        '''
        bodiesSerialized=[]
        for body in bodies:
            x = body.position.x
            y = body.position.y
            z = body.position.z
            
            vx=body.velocity.x
            vy=body.velocity.y
            vz=body.velocity.z
            
            pos = np.array([x,y,z])
            vel = np.array([vx,vy,vz])
            
            bodiesSerialized.append([vel,pos])

        f, time = self.solver.advance(bodiesSerialized,time)
            
#        print f.shape
#        print f
        for i in range(len(bodies)):
            
            bodies[i].velocity.x = f[0][i][0]
            bodies[i].velocity.y = f[0][i][1]
            bodies[i].velocity.z = f[0][i][2]
            
            bodies[i].position.x = f[1][i][0]
            bodies[i].position.y = f[1][i][1]
            bodies[i].position.z = f[1][i][2]
        
        return bodies, time
        
class nBody(object):
    '''Contains methods for use in calculating positions and velocities of
    bodies under the influence of a central gravity object.

    Attributes
    -----
    solver : object from Solver class
        Determines the type of approximation you want.
    M : float
        The mass of the central gravity object.
    '''
    def __init__(self,solver,M):
        self.G = -6.67408*10**(-11) # m^3/(kgs^2), the Gravitational Constant

        self.M = M
        self.solver = solver
        solver.diffEq = self.diffEq
    
    def diffEq(self,F,t):
        ''' Standard differential equations for motion.
        
        Parameters
        -----
        F : numpy array
            A numpy array containing two instances of a Vector class, first for
            velocity, second for position.

        Returns
        -----
        A numpy array containing the velocity and position Vectors as two lists
        within the array.
        '''
        G = self.G
        M = self.M
        velArray=F[0]
        posArray=F[1]
        accArray=np.empty([0,3])
        dpdtArray=np.empty([0,3])

#        print 'velArray = ' + repr(velArray)
#        print 'posArray = ' + repr(posArray)
            
        wholeLottaStuff=[]
            
        for j in range(len(F)):
            me=posArray[j,:]
            them=np.delete(posArray,j,axis=0)
            masses=np.delete(self.masses,j,axis=0)
            diff=them-me
            
#            print 'posArray = ' + repr(posArray)
#            print 'me = ' + repr(me)
#            print 'them = ' + repr(them)
#            print 'self.masses = ' + repr(self.masses)
#            print 'masses = ' + repr(masses)
#            print 'diff = ' +repr(diff)            
            
            rc=np.sum(diff**2,axis=1)**(-3/2)
            den=(np.ones(3*len(self.masses)).reshape(len(self.masses),3)*rc).T
            dvdt=G*masses*diff*den
        
            

#        print 'acc = ' + repr(dvdt)
#        print 'vel = ' + repr(velArray)
        return np.array([dvdt,velArray])
    
    def advance(self,bodies,time):
        '''Advances the given approximation forward by a step.
        
        Parameters
        -----
        body : an instance of GravBody
            The body whose positions and velocities are being calculated.
        time : float
            The time at which the approximation is being advanced.
            
        Returns
        -----
        body : an instance of GravBody
            The object whose positions and velocities are being calculated.
        time : float
            The new time after having advanced a step.
        '''
        
        posArray=np.empty([0,3])
        velArray=np.empty([0,3])
#        print 'velArray = ' + repr(velArray)
#        print 'posArray = ' + repr(posArray)
        self.masses=np.array([])
        for body in bodies:
            x = body.position.x
            y = body.position.y
            z = body.position.z
            m = body.mass
            
            
            vx=body.velocity.x
            vy=body.velocity.y
            vz=body.velocity.z
            
            pos = np.array([x,y,z])
            vel = np.array([vx,vy,vz])
            self.masses=np.append(self.masses,m)
            posArray = np.vstack([posArray,pos])
            velArray = np.vstack([velArray,vel])
            
#            print 'velArray = ' + repr(velArray)
#            print 'posArray = ' + repr(posArray)


        f, time = self.solver.advance(np.array([velArray,posArray]),time)
            
#        print f.shape
#        print f
        for i in range(len(bodies)):
            
            bodies[i].velocity.x = f[0][i][0]
            bodies[i].velocity.y = f[0][i][1]
            bodies[i].velocity.z = f[0][i][2]
            
            bodies[i].position.x = f[1][i][0]
            bodies[i].position.y = f[1][i][1]
            bodies[i].position.z = f[1][i][2]
        
        return bodies, time
        
        
def main():

    g = -9.81
    step = 0.1 
    
    vSim=Simulation.FindInitialV(Solver.RK2(step,UniformGravity.diffEq),g,150.,0)
    
    vSearch=Search.Bisect(vSim.vError,1000.,0.00001)
    actualX = vSearch.find()
    print 'Bisect returns ' + repr(actualX)
    
    vSearch=Search.Newton(vSim.vError,1000.,0.00001,0.001)
    actualX = vSearch.find()
    print "Newton's method returns " + repr(actualX)
    
if __name__ == '__main__': main()