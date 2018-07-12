# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 16:22:57 2016

@author: John
"""
import Body
import Physics
import Solver
import numpy as np
import copy as cp
import Vector
import math
import matplotlib.pyplot as plt

class CoolingSim(object):
    '''Use to simulate through a cooling curve approximation.
    
    Attributes
    -----
    solver : a solver class
        Used to determine the type of approximation.
    Ta : float
        The ambient temperature used for calculating the cooling curve.
    k : float
        The cooling constant.
    body : a ThermalBody object
        The body for which the cooling curve is being calculated.
    stopCondition : function
        A function that takes in a current y value (body.temperature), an 
        initial y value, and a final y value (usually Ta) and returns a boolean
        value depending on whether the desired conditions are met.
    '''
    def __init__(self,solver,Ta,k,body,stopCondition):
        
        self.solver = solver
        self.Ta = Ta
        self.k = k
        self.body = body
        self.stopCondition = stopCondition
        self.time = 0
        self.t = [self.time]
        self.T = [self.body.temperature]
        self.initialT = body.temperature
    
    def simulate(self):
        '''Uses a while loop to simulate through a cooling curve.
        
        Parameters
        -----
        None
        
        Returns
        -----
        t : list
            The list of times at which the temperatures were calculated.
        T : list
            The list of temperatures
        '''

        physics = Physics.NewtonCooling(self.solver,self.Ta,self.k)
        
        while self.stopCondition(self.body.temperature,self.initialT,self.Ta,) == False:
            self.body, self.time = physics.advance(self.body,self.time)
            self.t.append(self.time)
            self.T.append(self.body.temperature)
            
        return self.t, self.T
        
class TrajectorySim(object):
    '''Use to simulate the trajectory of a GravBody object.
    
    Attributes
    -----
    solver : an instance of a solver class
        Used to determine the type of approximation.
    bodies : list of instances of GravBody
        The bodies whose trajectories are being calculated.
    stopCondition : function
        A function that takes in a current velocity and returns a boolean
        value depending on whether the desired conditions are met (False until
        the condition is met and becomes True).
    M : float
        The mass of the central object. This is only used for inverse square or
        central body orbits.
    a : instance of Vector
        The acceleration vector due to gravity, defaults at x=0, y=0, and 
        z=-9.81. This is only used for a standard, simplified trajectory.
    c : float
        The coefficient due to drag. Leave at 0 for no drag.
    '''
    def __init__(self,solver,bodies,stopCondition,M,a=Vector.Vector(x=0,y=0,z=-9.81)
                                                    ,c=0):
        self.G = -6.67408*10**(-11) # m^3/(kgs^2), the Gravitational Constant
        self.solver = solver
        self.a = a
        self.M = M
        self.c = c

        self.bodies = bodies
        self.stopCondition = stopCondition
        self.time = 0
        self.t = [self.time]
        self.bodiesList=[]
        
        # Given a number of bodies at an initial time frame, I create a list
        # corresponding to each body into which I put all the future frames
        # of that body, and then put each of those lists into a new list that
        # holds them all together. This way, if I want all the information 
        # regarding the first body, I can just do bodiesList[0].
        for i in range(len(bodies)):
            subBodiesList = [cp.deepcopy(self.bodies[i])]
            self.bodiesList.append(subBodiesList)
            
        self.bodiesListInit = cp.deepcopy(self.bodiesList)
    
    def simulate(self):
        '''Uses a while loop to simulate through a trajectory.
        
        Parameters
        -----
        None
        
        Returns
        -----
        t : list
            The list of times at which the velocities and positions were calculated.
        f : list
            A list consisting of snapshots of the GravBody being simulated at
            every point in the trajectory.
        '''

        physics = Physics.UniformGravity(self.solver,self.a,c=self.c)
        
        while self.stopCondition(self.body.position.z) == False:
            self.body, self.time = physics.advance(self.body,self.time)
            self.t.append(self.time)
#            self.height.append(self.body.height)
#            self.velocity.append(self.body.velocity)
            self.bodyList.append(cp.deepcopy(self.body))

        return self.t, self.bodyList
        
    def simulateInvSqr(self):
        '''Uses a while loop to simulate through a trajectory that uses 
        gravity inversely proportional to radius squared.
        
        Parameters
        -----
        None
        
        Returns
        -----
        t : list
            The list of times at which the velocities and positions were calculated.
        f : numpy array
            An array consisting of a list of velocities and a list of heights
            in that order.
        '''

        physics = Physics.InverseSquareGravity(self.solver,self.M)
        
        while self.stopCondition(self.body.velocity) == False:
            self.body, self.time = physics.advance(self.body,self.time)
            self.t.append(self.time)
#            self.height.append(self.body.height)
#            self.velocity.append(self.body.velocity)
            self.bodyList.append(cp.deepcopy(self.body))


        return self.t, self.bodyList
        
class FindInitialAngle(object):
    '''Use to find the initial angle of launch for a GravBody object. Note that
    this requires that you have to have a stopCondition function in your
    Physics.py file.
    
    Attributes
    -----
    solver : an instance of a solver class
        Used to determine the type of approximation.
    a : instance of Vector
        The acceleration vector due to gravity.
    '''
    def __init__(self,solver,a=Vector.Vector(x=0,y=0,z=-9.81),c=0):
        self.solver = solver
        self.a = a
        self.c = c
        
    def angleError(self, phi):
        '''The error function used in a Search method for finding an angle that
        gives an equal max height and distance.
        
        Parameters
        -----
        phi : float
            The angle whose error is being calculated.
            
        Returns
        -----
        float
            The difference between the max height and the distance traveled for
            the given phi.
        '''
        rock = Body.GravBody(Vector.Vector(r=100.,theta=0,phi=phi),Vector.Vector(x=0,y=0,z=0))
        
        rockTrajectory = TrajectorySim(self.solver,rock,Physics.stopCondition,c=self.c)
        t,bodies = rockTrajectory.simulate()
        
        z=[b.position.z for b in bodies]
#        y=[b.position.y for b in bodies]
        x=[b.position.x for b in bodies]
#        print max(z)
#        print x[-1]
#        
#        plt.plot(x,z)
        
        dist=x[-1]
        
        return max(z)-dist
        
class FindInitialV(object):
    '''Use to find the velocity of a GravBody object. Note that this requires
    that you have to have a stopCondition function in your Physics.py file.
    
    Attributes
    -----
    solver : an instance of a solver class
        Used to determine the type of approximation.
    g : float
        The acceleration due to gravity (should be negative).
    desiredHeight: float
        The height you want to the velocity required to reach.
    startHeight : float
        The height your object starts at.
    '''
    
    def __init__(self,solver,g,desiredHeight,startHeight):
        self.solver = solver
        self.g = g
        self.desiredHeight = desiredHeight
        self.vGuess = 1000.
        self.hGuess = 51000.
        self.startHeight = startHeight
        
    def vError(self, v):
        '''The error function used in a Search method for finding an initial
        velocity that gives a trajectory that reaches the desired height.
        
        Parameters
        -----
        v : float
            The velocity whose error is being calculated.
            
        Returns
        -----
        float
            The difference between the calculated height and the desired height.
        '''
        step=0.01
        
        skull = Body.GravBody(Vector.Vector(x=0,y=0,z=v),Vector.Vector(x=0,y=0,z=self.startHeight))
        solver = Solver.RK2(step,Physics.UniformGravity.diffEq)
        
        skullFlight = TrajectorySim(solver,Vector.Vector(x=0,y=0,z=self.g),skull,Physics.stopCondition)
        t,bodies = skullFlight.simulate()
        
        hGuess = bodies[-1].position.z
        
        return hGuess-self.desiredHeight  
            
    def findV(self,tolerance):
        '''Finds the velocity given a desired height within a given tolerance
        for a body on Earth.
        
        Parameters
        -----
        tolerance : float
            How close you want to be to the theoretically ideal value.
            
        Returns
        -----
        The estimated velocity.
        '''
        
        desiredHeight=self.desiredHeight
        vGuess=self.vGuess
        hGuess=self.hGuess
        g = -9.81
        step = 0.1 
        lastvGuess = vGuess
        up = False
        down = False
        
        while abs(hGuess-desiredHeight)>tolerance:
            
            skull = Body.GravBody(vGuess,self.startHeight)
            solver = Solver.RK2(step,Physics.UniformGravity.diffEq)
            
            skullFlight = TrajectorySim(solver,g,skull,Physics.stopCondition)
            t,vAndh = skullFlight.simulate()
            
            hGuess = vAndh[1][-1]
            
            if hGuess-desiredHeight < 0:
                up = True
                if down == False:
                    lastvGuess = vGuess
                    vGuess = vGuess*2
                else:
                    lvgHold = vGuess
                    vGuess = vGuess + abs(lastvGuess-vGuess)*0.5
                    lastvGuess = lvgHold
                    
            elif hGuess-desiredHeight > 0:
                down = True
                if up == False:
                    lastvGuess = vGuess
                    vGuess = vGuess/2.
                else:
                    lvgHold = vGuess
                    vGuess = vGuess - abs(vGuess-lastvGuess)*0.5
                    lastvGuess = lvgHold
            else: print('Error?')
        
        return vGuess
        
class OrbitSim(TrajectorySim):
    '''Use to simulate the orbit of a GravBody object.
    
    Attributes
    -----
    solver : an instance of a solver class
        Used to determine the type of approximation.
    bodies : list of instances of GravBody
        The bodies whose trajectories are being calculated.
    stopCondition : function
        A function that takes in a current velocity and returns a boolean
        value depending on whether the desired conditions are met (False until
        the condition is met and becomes True).
    M : float
        The mass of the central object. This is only used for inverse square or
        central body orbits.
    orbits : int
        The desired number of orbits through which you want to simulate.
    '''
    def __init__(self,solver,bodies,stopCondition,M,physics):
        super(OrbitSim,self).__init__(solver,bodies,stopCondition,M)
        self.physics=physics
        self.desiredOrbitRadius=None
        self.numberOrbits=0
        
    def simulateOrbit(self):
        '''Uses a while loop to simulate through an orbit that uses 
        gravity inversely proportional to radius squared.
        
        Parameters
        -----
        None
        
        Returns
        -----
        t : list
            The list of times at which the velocities and positions were calculated.
        f : numpy array
            An array consisting of a list of velocities and a list of heights
            in that order.
        '''

        physics = self.physics(self.solver,self.M)

        while self.stopCondition(self) == False:
            self.bodies, self.time = physics.advance(self.bodies,self.time)
            self.numberOrbits=self.countOrbits(self.bodiesList,self.numberOrbits)
            self.t.append(self.time)
            for i in range(len(self.bodies)):
                self.bodiesList[i].append(cp.deepcopy(self.bodies[i]))        

        return self.t, self.bodiesList       
        
    def period(self,time,bodies):
        '''Calculates the orbital period of a body.
        
        Parameters
        -----
        time : list
            List of times corresponding to the snapshots taken for each of the
            bodies.
        bodies : list of GravBody objects
            List of GravBody objects corresponding to each snapshot in an
            orbit.
        '''
#        return 2.*math.pi/math.sqrt(self.G*self.M/self.body.position.r**3)
    
        period=0
        
        if self.orbits==1:
            period = time[-1]-time[0]
        else:

            originalPos = bodies[0].position
            stopSpan=abs(bodies[1].position-originalPos)
            
            for n in range(3,len(bodies)):
                
                currentPos = bodies[n].position
                diff = abs(currentPos-originalPos)
                
                if diff.r<=stopSpan.r/2.:
                    period = time[n]-time[0]
                    break
                
        return period
        
    def semiMajorAxis(self,bodies):
        '''Calculates the semimajor axis of a body's orbit.
        
        Parameters
        -----
        bodies : list of GravBody objects
            List of GravBody objects corresponding to each snapshot in an
            orbit.
        '''
#        return 1./((2./self.body.position.r)-(self.body.velocity.r**2/(self.G*self.M)))


        r=[g.position.r for g in self.bodies[0]]        

        rmax=max(r)
        rmin=min(r)
        
        return (rmax+rmin)/2.
        
    def totalEnergy(self,body):
        PE = abs(self.G*self.M*body.mass/body.position.r)
        KE = abs(body.velocity.r**2*0.5*body.mass)
        return PE+KE
        
    @staticmethod   
    def countOrbits(bodiesList,numberOrbits):
        '''This takes in orbital info and a desired number of orbits and will 
        count the number of orbits achieved.
        
        Parameters
        -----
        bodies : array of GravBody objects
            All of the bodies at each timestep in a simulation.
        desiredOrbits : int
            The number of orbits to simulate through in total.
        numOrbits : int
            The number of orbits that the simulation has gone through thus far.
        '''

        bodies=bodiesList[0]
        if len(bodies)<3:
            return numberOrbits
#        elif len(bodies)>50000:
#            print "Didn't reach desired number of orbits!"
#            return desiredOrbits
        else:
            originalPos = bodies[0].position
            stopSpan=abs(bodies[1].position-originalPos)
            currentPos = bodies[-1].position
            diff = abs(currentPos-originalPos)
            if diff.r <= stopSpan.r/2.:
                numberOrbits=numberOrbits+1
    #            print numOrbits
#            phiShip=[g.position.theta for g in bodies[3:]]
#            phiShip=np.array(phiShip)
#            phiShip=phiShip-originalPos
#            numOrbits=sum(phi < stopSpan/2. for phi in phiShip)
            #Count number of phiShip<stopSpan
            return numberOrbits
        
    def vError(self, v):
        '''The error function used in a Search method for finding an initial
        velocity that gives a trajectory that reaches the desired height.
        
        Parameters
        -----
        v : float
            The velocity whose error is being calculated.
            
        Returns
        -----
        float
            The difference between the calculated height and the desired height.
        '''
        print('Velocity being tested ' + repr(v))
        
        self.time = 0
        self.t = [self.time]

        self.bodiesList=[]
        self.bodies[0].velocity.y=v
        for i in range(len(self.bodies)):
            subBodiesList = [cp.deepcopy(self.bodies[i])]
            self.bodiesList.append(subBodiesList)

        t,bodies = self.simulateOrbit()
        
        rShip=[g.position.r for g in bodies[0]]  
        
        yShip=[g.position.y for g in bodies[0]]
        xShip=[g.position.x for g in bodies[0]]

        
        plt.plot(xShip,yShip,'.')

#        plt.plot(0,0, 'y*')
#        plt.xlabel('Distance (m)')
#        plt.ylabel('Distance (m)')
#        plt.title('Orbits Yay!')
#        
        rGuess=max(rShip)
        vError=rGuess-self.desiredRadius
#        print 'vError = ' + repr(vError)
        return vError
        
    def tError(self, t):
        '''The error function used in a Search method for finding an initial
        velocity that gives a trajectory that reaches the desired height.
        
        Parameters
        -----
        v : float
            The velocity whose error is being calculated.
            
        Returns
        -----
        float
            The difference between the calculated height and the desired height.
        '''
        print('Time being tested ' + repr(t))

        self.time = 0
        self.t = [self.time]
        self.bodiesList=[]
        for i in range(len(self.bodies)):
            subBodiesList = [cp.deepcopy(self.bodies[i])]
            self.bodiesList.append(subBodiesList)

        refTime,refBodies = self.simulateOrbit()
        
        newStart = refTime.index(t)
        self.bodiesList = refBodies[newStart]
        
        t,bodies=self.simulateOrbit()

        posError=min(bodies[0].position-bodies[2].position)
#        print 'vError = ' + repr(vError)
        return posError.position.r
        
class BinarySim(TrajectorySim):
    '''Use to simulate the orbit of a binary system of GravBody objects.
    
    Attributes
    -----
    solver : an instance of a solver class
        Used to determine the type of approximation.
    bodies : list of instances of GravBody
        The bodies whose trajectories are being calculated.
    stopCondition : function
        A function that takes in a current velocity and returns a boolean
        value depending on whether the desired conditions are met (False until
        the condition is met and becomes True).
    physics : instance of a physics class
        Used to determine the diffEQ that solver uses.
    M1 : float
        Mass of body 1.
    M2 : float
        Mass of body 2.
    a1 : float
        Semimajor axis of body 1.
    e : float between 0 and 1
        The eccentricity of the bodies' orbits.
    '''
    def __init__(self,solver,stopCondition,physics,M1,M2,a1,e):

        self.M1=M1
        self.M2=M2
        self.a1=a1
        self.e=e
        self.physics=physics
        self.desiredOrbitRadius=None
        self.numberOrbits=0
        self.G=6.67408*10**-11
        
        self.v1=math.sqrt(self.G*M2**3*(1+e)/(a1*(M1+M2)**2*(1-e)))
        self.v2=-(M1/M2)*self.v1
        self.r1=a1-a1*e
        self.r2=-(M1/M2)*self.r1
        
        M=1.989*10**30
        
        body1x=self.r1
        body1y=0
        body1z=0
        body1pos=Vector.Vector(x=body1x,y=body1y,z=body1z)
        
        body2x=self.r2
        body2y=0
        body2z=0
        body2pos=Vector.Vector(x=body2x,y=body2y,z=body2z)
        
        body1vx=0
        body1vy=self.v1
        body1vz=0
        body1vel=Vector.Vector(x=body1vx,y=body1vy,z=body1vz)
        
        body2vx=0
        body2vy=self.v2
        body2vz=0
        body2vel=Vector.Vector(x=body2vx,y=body2vy,z=body2vz)
        
        self.body1=Body.GravBody(body1vel,body1pos,M1)
        self.body2=Body.GravBody(body2vel,body2pos,M2)
        bodies=[self.body1,self.body2]
        super(BinarySim,self).__init__(solver,bodies,stopCondition,M)

    def simulateOrbit(self):
        '''Uses a while loop to simulate through an orbit that uses 
        gravity inversely proportional to radius squared.
        
        Parameters
        -----
        None
        
        Returns
        -----
        t : list
            The list of times at which the velocities and positions were calculated.
        f : numpy array
            An array consisting of a list of velocities and a list of heights
            in that order.
        '''

        physics = self.physics(self.solver,self.M)

        while self.stopCondition(self) == False:
            self.bodies, self.time = physics.advance(self.bodies,self.time)
            self.numberOrbits=self.countOrbits(self.bodiesList,self.numberOrbits)
            self.t.append(self.time)
            for i in range(len(self.bodies)):
                self.bodiesList[i].append(cp.deepcopy(self.bodies[i]))        

        return self.t, self.bodiesList     
            
    @staticmethod   
    def countOrbits(bodiesList,numberOrbits):
        '''This takes in orbital info and a desired number of orbits and will 
        count the number of orbits achieved.
        
        Parameters
        -----
        bodies : array of GravBody objects
            All of the bodies at each timestep in a simulation.
        desiredOrbits : int
            The number of orbits to simulate through in total.
        numOrbits : int
            The number of orbits that the simulation has gone through thus far.
        '''

        bodies=bodiesList[0]
        if len(bodies)<3:
            return numberOrbits
#        elif len(bodies)>50000:
#            print "Didn't reach desired number of orbits!"
#            return desiredOrbits
        else:
            originalPos = bodies[0].position
            stopSpan=abs(bodies[1].position-originalPos)
            currentPos = bodies[-1].position
            diff = abs(currentPos-originalPos)
            if diff.r <= stopSpan.r/2.:
                numberOrbits=numberOrbits+1
    #            print numOrbits
#            phiShip=[g.position.theta for g in bodies[3:]]
#            phiShip=np.array(phiShip)
#            phiShip=phiShip-originalPos
#            numOrbits=sum(phi < stopSpan/2. for phi in phiShip)
            #Count number of phiShip<stopSpan
            return numberOrbits
            
class ExoSim(TrajectorySim):
    '''Use to simulate the orbit of a binary system of GravBody objects.
    
    Attributes
    -----
    solver : an instance of a solver class
        Used to determine the type of approximation.
    bodies : list of instances of GravBody
        The bodies whose trajectories are being calculated.
    stopCondition : function
        A function that takes in a current velocity and returns a boolean
        value depending on whether the desired conditions are met (False until
        the condition is met and becomes True).
    physics : instance of a physics class
        Used to determine the diffEQ that solver uses.
    Ms : float
        Mass of star.
    Mp : float
        Mass of planet.
    Rs : float
        The radius of the star.
    Rp : float
        The radius of the planet.
    ap : float
        Semimajor axis of the planet's orbit.
    e : float between 0 and 1
        The eccentricity of the bodies' orbits.
    omega : float
        The angle of periastron.
    i : float
        The inclination of the orbit.
    '''
    def __init__(self,solver,stopCondition,physics,Ms,Mp,Rs,Rp,ap,e,omega,i):
        
        self.Rs=Rs
        self.Rp=Rp
        self.omega=omega
        self.i=i
        self.Ms=Ms
        self.Mp=Mp
        self.ap=ap
        self.e=e
        self.physics=physics
        self.desiredOrbitRadius=None
        self.numberOrbits=0
        self.G=6.67408*10**-11
        
        self.vp=math.sqrt(self.G*Ms**3*(1+e)/(ap*(Mp+Ms)**2*(1-e)))
        self.vs=-(Mp/Ms)*self.vp
        self.rp=ap-ap*e
        self.rs=-(Mp/Ms)*self.rp
        
        M=1.989*10**30
        
        body1x=self.rp
        body1y=0
        body1z=0
        body1pos=Vector.Vector(x=body1x,y=body1y,z=body1z)
        body1pos.rotZ(omega)
        body1pos.rotX(i)
        
        body2x=self.rs
        body2y=0
        body2z=0
        body2pos=Vector.Vector(x=body2x,y=body2y,z=body2z)
        body2pos.rotZ(omega)
        body2pos.rotX(i)
        
        body1vx=0
        body1vy=self.vp
        body1vz=0
        body1vel=Vector.Vector(x=body1vx,y=body1vy,z=body1vz)
        body1vel.rotZ(omega)
        body1vel.rotX(i)
        
        body2vx=0
        body2vy=self.vs
        body2vz=0
        body2vel=Vector.Vector(x=body2vx,y=body2vy,z=body2vz)
        body2vel.rotZ(omega)
        body2vel.rotX(i)
        
        self.body1=Body.GravBody(body1vel,body1pos,Mp)
        self.body2=Body.GravBody(body2vel,body2pos,Ms)
        bodies=[self.body1,self.body2]
        
        super(ExoSim,self).__init__(solver,bodies,stopCondition,M)

    def simulateOrbit(self):
        '''Uses a while loop to simulate through an orbit that uses 
        gravity inversely proportional to radius squared.
        
        Parameters
        -----
        None
        
        Returns
        -----
        t : list
            The list of times at which the velocities and positions were calculated.
        f : numpy array
            An array consisting of a list of velocities and a list of heights
            in that order.
        '''
        self.AeclipsedList=[]
        physics = self.physics(self.solver,self.M)
        Rs=self.Rs
        Rp=self.Rp

        while self.stopCondition(self) == False:
            self.bodies, self.time = physics.advance(self.bodies,self.time)
            
            d=np.sqrt((self.bodies[0].position.x-self.bodies[1].position.x)**2
                +(self.bodies[0].position.y-self.bodies[1].position.y)**2)
            if d>=(self.Rp+self.Rs):
                Aeclipsed=0
            elif self.bodies[0].position.z>0:
                Aeclipsed=0
            elif d<=self.Rs-self.Rp:
                Aeclipsed=np.pi*self.Rp**2
            else:
                x=(Rp**2+d**2-Rs**2)/(2*d)
                h=np.sqrt(Rp**2-x**2)
                ApTri=  x  *h
                AsTri=(d-x)*h
                thetap=np.arccos(x/Rp)
                thetas=np.arccos((d-x)/Rs)
                ApSec=thetap*Rp**2
                AsSec=thetas*Rs**2
                Ap=ApSec-ApTri
                As=AsSec-AsTri
                Aeclipsed=Ap+As
            self.AeclipsedList.append(Aeclipsed)
            
            self.numberOrbits=self.countOrbits(self.bodiesList,self.numberOrbits)
            self.t.append(self.time)
            for i in range(len(self.bodies)):
                self.bodiesList[i].append(cp.deepcopy(self.bodies[i]))        

        return self.t, self.bodiesList     
            
    @staticmethod   
    def countOrbits(bodiesList,numberOrbits):
        '''This takes in orbital info and a desired number of orbits and will 
        count the number of orbits achieved.
        
        Parameters
        -----
        bodies : array of GravBody objects
            All of the bodies at each timestep in a simulation.
        desiredOrbits : int
            The number of orbits to simulate through in total.
        numOrbits : int
            The number of orbits that the simulation has gone through thus far.
        '''

        bodies=bodiesList[0]
        if len(bodies)<3:
            return numberOrbits
#        elif len(bodies)>50000:
#            print "Didn't reach desired number of orbits!"
#            return desiredOrbits
        else:
            originalPos = bodies[0].position
            stopSpan=abs(bodies[1].position-originalPos)
            currentPos = bodies[-1].position
            diff = abs(currentPos-originalPos)
            if diff.r <= stopSpan.r/2.:
                numberOrbits=numberOrbits+1
    #            print numOrbits
#            phiShip=[g.position.theta for g in bodies[3:]]
#            phiShip=np.array(phiShip)
#            phiShip=phiShip-originalPos
#            numOrbits=sum(phi < stopSpan/2. for phi in phiShip)
            #Count number of phiShip<stopSpan
            return numberOrbits
            
    def lightCurve(self):
        '''For use in plotting a light curve.
        
        Parameters
        -----
        none
        
        Returns
        -----
        An array of ratios that can be used to plot a light curve.
        '''
        As = np.pi*self.Rs**2
        return (As-np.array(self.AeclipsedList))/As
        
    @staticmethod
    def lightCurveStatic(Rs,AeclipsedList):
        '''For use plotting a light curve outside of a simulation.
        '''
        As = np.pi*Rs**2
        return (As-np.array(AeclipsedList))/As
        
    @staticmethod
    def AEclipsed(dArray,Rp,Rs ):
        '''For use plotting a light curve outside of a simulation.
        '''
        AeclipsedArray=[]
        for d in dArray:
            if d>=(Rp+Rs):
                Aeclipsed=0
#            elif sim.bodies[0].position.z>0:
#                Aeclipsed=0
            elif d<=Rs-Rp:
                Aeclipsed=np.pi*Rp**2
            else:
                x=(Rp**2+d**2-Rs**2)/(2*d)
                h=np.sqrt(Rp**2-x**2)
                ApTri=  x  *h
                AsTri=(d-x)*h
                thetap=np.arccos(x/Rp)
                thetas=np.arccos((d-x)/Rs)
                ApSec=thetap*Rp**2
                AsSec=thetas*Rs**2
                Ap=ApSec-ApTri
                As=AsSec-AsTri
                Aeclipsed=Ap+As
                print('thetap = ' + repr(thetap))
                print('thetas = ' + repr(thetas))


            AeclipsedArray.append(Aeclipsed)
            
        curve = ExoSim.lightCurveStatic(Rs,AeclipsedArray)
        plt.plot(curve,'.')
        return AeclipsedArray
        
class SolarModel(object):
    '''
    '''
    def __init__(self, stepsize, solver, physics, stopCondition):

        self.G = 6.67408*10**(-11) # m^3/(kgs^2), the Gravitational Constant
        self.solver = solver
        self.physics = physics

        self.stopCondition = stopCondition
        self.time = 0
        self.t = [self.time]
        self.originalStepSize=cp.deepcopy(self.solver.stepsize)

#        self.bodiesList=[]
        
    def simulate(self):
        '''
        '''
        physics = self.physics

        while self.stopCondition(self) == False:
            self.bodies, self.time = physics.advance(self.bodies,self.time)
            self.solver.stepsize=self.originalStepSize

#            self.numberOrbits=self.countOrbits(self.bodiesList,self.numberOrbits)
            self.t.append(self.time)     

        return self.t, self.bodies  
        
    @staticmethod
    def closed2BodyAutoBuilder(omega,i,mBig,mSmall,ap,e):
        '''
        '''
        G = 6.67408*10**(-11)
        vp=math.sqrt(G*mBig**3*(1+e)/(ap*(mSmall+mBig)**2*(1-e)))
        vs=-(mSmall/mBig)*vp
        rp=ap-ap*e
        rs=-(mSmall/mBig)*rp
        
        body1x=rp
        body1y=0
        body1z=0
        body1pos=Vector.Vector(x=body1x,y=body1y,z=body1z)
        body1pos.rotZ(omega)
        body1pos.rotX(i)
        
        body2x=rs
        body2y=0
        body2z=0
        body2pos=Vector.Vector(x=body2x,y=body2y,z=body2z)
        body2pos.rotZ(omega)
        body2pos.rotX(i)
        
        body1vx=0
        body1vy=vp
        body1vz=0
        body1vel=Vector.Vector(x=body1vx,y=body1vy,z=body1vz)
        body1vel.rotZ(omega)
        body1vel.rotX(i)
        
        body2vx=0
        body2vy=vs
        body2vz=0
        body2vel=Vector.Vector(x=body2vx,y=body2vy,z=body2vz)
        body2vel.rotZ(omega)
        body2vel.rotX(i)
        
        body1=Body.GravBody(body1vel,body1pos,mSmall)
        body2=Body.GravBody(body2vel,body2pos,mBig)
        bodies=[body1,body2]
        
        return bodies
        
    @staticmethod
    def manualOrbitBuilder(xPos,yPos,zPos,xVel,yVel,zVel,mass):
        '''
        '''
        
        bodies=[]
        
        for i in range(len(xPos)):
            bodyPos = Vector.Vector(x=xPos[i],y=yPos[i],z=zPos[i])
            bodyVel = Vector.Vector(x=xVel[i],y=yVel[i],z=zVel[i])
            body    = Body.GravBody(bodyVel,bodyPos,mass[i])
            bodies.append(body)
            
        return bodies
            