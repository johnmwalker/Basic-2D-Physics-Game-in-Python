"""
Created on Tue Mar 08 15:50:22 2016

This is the final draft of the game I made for PHYS 398 in Spring 2016
(sophomore year). I have not bothered to update it much since then. See other 
projects for better coding.

Goal of game: 
    Suck the energy out of the planets to make them fall into the black hole.
    You'll lose if you fall in yourself.

Controls: 
    Numpad + and - speed up and slow down the simulation
    > and < (. and ,) are "cheats" to give your thruster more power.
    Spacebar activates the "energy sucking" function. Planets in the crosshairs
        will slow down.
    Arrow keys control motion of ship.
    
Issues:
    I couldn't figure out how to fix the rotation animation of the ship, so it
    doesn't really rotate about the center of the image. I think this may have
    been due to a bug in PyGame.
    
    Collisions were poorly implemented. For some reason I decided to use trig
    to calculate forces. Later projects built on this code use vectors for
    greatly improved accuracy.
    
    This was my first big project, so many coding conventions were likely
    ignored.
"""
import Physics
import Simulation
import Solver
import numpy as np
import StopConditions
import pygame as pg

# Define some colors
BLACK    = (   0,   0,   0)
WHITE    = ( 255, 255, 255)
GREEN    = (   0, 255,   0)
RED      = ( 255,   0,   0)

# Choose the "level" (defines starting conditions)
# Level 1: Standard game layout. Several planets (one of them has a tiny moon!)
# Level 2: Roughly to-scale-ish solar system
# Level 3: Planet sandbox to demonstrate collisions
level=1

class App(object):
    ''' Core class that connects all the dots.
    '''
    def __init__(self, level):
        
        if level==1:
            self.scale=1.*10**-8.8 # -9 for level1
            self.sPT=10. # 10 for level1
            pg.init()       # Initialize pygame 
            self.initializeSimulation()
            self.level1()
        elif level==2:
            self.scale=1.*10**-8.2 # -8.3 or less for level2
            self.sPT=50. # 100 for level2
            pg.init()       # Initialize pygame 
            self.initializeSimulation()
            self.level2()
        elif level==3:
            self.scale=1.*10**-6 # -6 for level3
            self.sPT=1000. # 1000 for level3
            pg.init()       # Initialize pygame 
            self.initializeSimulation()
            self.level3()
        
        ##### Finalize #####
        self.sprites = pg.sprite.RenderUpdates(self.listOfSprites) 
        self.sprites.clear(self.screen,self.background)
        self.thingThatGetsDrawn = self.sprites.draw(self.screen)
        pg.display.update(self.thingThatGetsDrawn)
        
    def level1(self):
        '''
        '''
        self.level=1
        ##### Parameters for Autobuilder #####

        M         =2.*10**34
        smallM    =[2.*10**27, 7.*10**22]
        a         =[4.0*10**11, 3.0*10**10.3]
        e         =[0,0]
        imageString = ["ganymedeSmall.png","gasPlanetSmall.png"]
        
        ##### Parameters for Manual builder #####

        mass=[M, 1*10**23  , 1*10**24 , 1*10**25 , 1*10**26 ]
        xPos=[0, 1*10**11, -1*10**11, 2.*10**11, 3.*10**11]
        yPos=[0, 0,          0,         0,         0        ]
        zPos=[0, 0,          0,         0,         0        ]
        
        vel=self.orbVel(M,np.array(mass[1:]),np.array(xPos[1:]))
        
        xVel=[0, 0,          0,         0,         0        ]
        yVel=[0, vel[0],     -vel[1],    vel[2],    vel[3]   ]
        zVel=[0, 0,          0,         0,         0        ]
        imageString2 =["blackhole.png","planetSmall.png","planet2Small.png",
                                       "planet3Small.png","planet4Small.png"]
        
        ##### Build Rocket #####
        
        mRocket = [6.*10**30]
        xRocket,yRocket,zRocket    = [5.0*10**11], [0], [0]
        vxRocket,vyRocket,vzRocket = [0],[self.orbVel(M,mRocket[0],
                                      xRocket[0])],[0]
        rocketImage = ["shipSmall.png"]
        self.rocketImageOrig=pg.image.load(rocketImage[0])
        
        self.genSpriteAndBodManual(xRocket,yRocket,zRocket,
                                   vxRocket,vyRocket,vzRocket,
                                   mRocket,rocketImage)
        self.rocketOnImage=pg.image.load("shipSmallThrusting.png")
        self.rocketOnImageOrig=pg.image.load("shipSmallThrusting.png")
                                   
        ##### Initialize ship characteristics #####
        self.listOfSprites[0].angle = 0
        self.listOfSprites[0].thrust = 20000.
        
        ##### Build bodies #####
        self.genSpriteAndBodManual(xPos, yPos, zPos, xVel, yVel, zVel, mass, 
                                                               imageString2)  
        self.genSpriteAndBod(smallM,e,M,a,imageString,
                        moon=True,keepCentral=False)
                        
        ##### Make central body a black hole and ship the 'main body' #####
        self.listOfSprites[1].blackHole = True
        self.listOfSprites[0].mainBody = True
        self.makeDetectionBox(self.whiteCrosshairs)
        
    def level2(self):
        ''' Solar System!
        '''
        self.level=2
        ##### Parameters for Manual builder #####
        mass=[1.989*10**30, 3.285*10**23, 4.867*10**24, 5.972*10**24, 6.39*10**23, 1.898*10**27, 5.683*10**26, 8.681*10**25, 1.024*10**26]
        xPos=[0,            57.9*10**9,   108.2*10**9,  149.6*10**9,  227.9*10**9, 778.3*10**9,  1427.0*10**9, 2871.0*10**9, 4497.1*10**9]
        yPos=[0,            0,            0,            0,            0,           0,            0,            0,            0,          ]
        zPos=[0,            0,            0,            0,            0,           0,            0,            0,            0,          ]
        
        vel=self.orbVel(1.989*10**30,np.array(mass[1:]),np.array(xPos[1:]))
#        print vel
        
        xVel=[0,            0,            0,            0,            0,           0,            0,            0,            0,          ]
        yVel=[0,            vel[0],       vel[1],       vel[2],       vel[3],      vel[4],       vel[5],       vel[6],       vel[7]      ]
        zVel=[0,            0,            0,            0,            0,           0,            0,            0,            0,          ]
        imageString =['sun.png','mercurySmall.png','venusSmall.png','earthSmall.png','marsSmall.png','jupiterSmall.png','saturn.png','uranusSmall.png','neptuneSmall.png']
        
        ##### Build Rocket #####
        mRocket = [2030.*10**3]
        xRocket,yRocket,zRocket    = [5.0*10**10], [0], [0]
        vxRocket,vyRocket,vzRocket = [0],[self.orbVel(1.989*10**30,mRocket[0],
                                      xRocket[0])],[0]
        rocketImage = ["shipSmall.png"]
        self.rocketImageOrig=pg.image.load(rocketImage[0])
        
        self.genSpriteAndBodManual(xRocket,yRocket,zRocket,
                                   vxRocket,vyRocket,vzRocket,
                                   mRocket,rocketImage)
        self.rocketOnImage=pg.image.load("shipSmallThrusting.png")
        self.rocketOnImageOrig=pg.image.load("shipSmallThrusting.png")
                                   
        ##### Initialize ship characteristics #####
        self.listOfSprites[0].angle= 0
        self.listOfSprites[0].thrust = 20000.
        
        ##### Build bodies #####
        self.genSpriteAndBodManual(xPos, yPos, zPos, xVel, yVel, zVel, mass, 
                                                               imageString)  
                        
        ##### Make central body a black hole and ship the 'main body' #####
        self.listOfSprites[1].blackHole = True
        self.listOfSprites[0].mainBody = True
        self.makeDetectionBox(self.whiteCrosshairs)

    def level3(self):
        '''
        '''
        self.level=3
        mass=[1.*10**22,1.*10**20,1.*10**19, 1.*10**22]
        xPos=[-5.*10**8,-2.*10**8,2.*10**8,5.*10**8]
        yPos=[0,0,2.*10**8,0]
        zPos=[0,0,0,0]
        xVel=[0,10.,0,0]
        yVel=[15.,10.,10.,-15.]
        zVel=[0,0,0,0]
        
        imageString = ['jupiterSmall.png','marsSmall.png','venusSmall.png','neptuneSmall.png']
        self.genSpriteAndBodManual(xPos,yPos,zPos,xVel,yVel,zVel,mass,imageString)
        self.listOfSprites[0].mainBody=True
        
    def advance(self):
        '''Advance the simulation one frame'''
            
        ship = self.listOfSprites[0]
        # Calculate step
        current_ticks = pg.time.get_ticks()
        self.step = (current_ticks - self.ticks)*self.secondsPerTick
        self.ticks = current_ticks

        # Run the simulation forward
        self.model.stopCondition=StopConditions.afterTimeSetup(time=self.t[-1]+
                                                                    self.step)
        self.t, self.bodies = self.model.simulate()  
        
#        print self.model.bodies[0].velocity
     
        # Update the display
#        coords = self.listOfSprites[0].rect
#        print len(self.t)
        if len(self.t) == 2: self.calcPosCoM()
#        print coords
        
        bgSize=[self.background.get_width(), self.background.get_height()]
        self.screen.blit(self.background,(self.BGcoords[0]-bgSize[0]/2.,
                                          self.BGcoords[1]-bgSize[1]/2.))
        self.sprites.update([ship.rect[0],ship.rect[1]],
            [ship.body.position.x,ship.body.position.y], [ship.image.get_width(), ship.image.get_height()],
            ship.angle)
        thingThatGetsDrawn = self.sprites.draw(self.screen)
        pg.display.update()
        
    def calcPosCoM(self):
        '''
        '''
        masses=np.array([s.body.mass for s in self.listOfSprites if s.body != None])
        xCoords  =np.array([s.rect[0] for s in self.listOfSprites if s.body != None])
        yCoords  =np.array([s.rect[1] for s in self.listOfSprites if s.body != None])
        
        x=np.sum(masses*xCoords)/np.sum(masses)
        y=np.sum(masses*yCoords)/np.sum(masses)
        
        self.BGcoords=[x,y]
        print(self.BGcoords)
        
        xVal  =np.array([s.body.position.x for s in self.listOfSprites if s.body != None])
        yVal  =np.array([s.body.position.y for s in self.listOfSprites if s.body != None])
        
        x=np.sum(masses*xVal)/np.sum(masses)
        y=np.sum(masses*yVal)/np.sum(masses)
        
        self.BGpos=[x,y]
        
        print('CoM coords calculated')
        
    def checkCollisions(self):
        ''' Checks for collisions between the necessary sprites.
        '''
        # Initialize type of collision
        circleCollide = pg.sprite.collide_circle_ratio(np.sqrt(2.8)/2.)
        
        # Check to see if we need to ignore the detection box
        if self.detectionBox==True:
            p=0
            q=1
        else: 
            p=1
            q=0
        
        # Build sprite group for checking collisions
        for i in range(len(self.listOfSprites)-q):
            if self.detectionBox==True: 
                spriteGroup = pg.sprite.Group(self.listOfSprites[i+1:-1])
            elif self.detectionBox==False: 
                spriteGroup = pg.sprite.Group(self.listOfSprites[i:])
                
            # Check collisions
            col = pg.sprite.spritecollide(self.listOfSprites[i],spriteGroup,
                                          dokill=False,collided=circleCollide)

            # Check if either colliding bodies are black holes. If so, delete
            # the necessary body.
            for j in range(i+1,len(self.listOfSprites)):
#                print j
                if self.listOfSprites[j] in col:
                    print('Body ' + repr(i) + ' and body ' + repr(j) +
                                                            ' collided.')
                    if self.listOfSprites[i].blackHole==True:
                        m=self.model.bodies[j].mass
                        self.model.bodies[i].mass += m # adds mass to black hole
                        self.listOfSprites[j].kill()
                        del self.listOfSprites[j]
                        del self.model.bodies[j]
                        if j==0:
                            self.gameOver=True
                        break
                    
                    elif self.listOfSprites[j].blackHole==True:
                        m=self.model.bodies[i].mass
                        self.model.bodies[j].mass += m
                        self.listOfSprites[i].kill()
                        del self.listOfSprites[i]
                        del self.model.bodies[i]
                        if i==0:
                            self.gameOver=True
                        break   
                    else: self.collision(self.model.bodies[i],
                                         self.model.bodies[j])
    
    def checkFloorCollisions(self):
        '''
        '''
        floorHeight=-5.*10**-9
        bodies=[s.body for s in self.listOfSprites if s.body!=None]
        for body in bodies:
            if body.position.y-floorHeight <= abs(body.velocity.y*self.secondsPerTick):
                body.velocity.y*=-1
                body.position.y= abs(body.velocity.y*self.secondsPerTick)-(body.position.y-floorHeight)
                print('Bounced off floor!')
                
                    
    def checkWin(self):
        ''' Checks to see if only you and other black holes are still around
        and triggers a win if so.
        '''
        # Check for detection box so we know to ignore it if necessary
        if self.detectionBox==True: p=1
        else: p=0
        
        # Count the number of mainbodies and black holes
        trueList=[]
        for sprite in self.listOfSprites:
            if sprite.mainBody == True or sprite.blackHole==True:
                trueList.append(True)
            else: break
        
        # If the number of black holes and main bodies is equal to the number
        # of bodies still around, you win!
        if len(trueList)==len(self.listOfSprites)-p:
            self.win = True
                   
    def collision(self,body1,body2):
        ''' Colliding spheres.
        '''
        
        theta1=body1.velocity.theta
        theta2=body2.velocity.theta
        rho=(body1.position-body2.position).theta
        m1=body1.mass
        m2=body2.mass
        v1=body1.velocity.r
        v2=body2.velocity.r
        
        v1x=((v1*np.cos(theta1-rho)*(m1-m2)+2*m2*v2*np.cos(theta2-rho))/
            (m1+m2)*np.cos(rho)+v1*np.sin(theta1-rho)*np.cos(rho+np.pi/2))
        v1y=((v1*np.cos(theta1-rho)*(m1-m2)+2*m2*v2*np.cos(theta2-rho))/
            (m1+m2)*np.sin(rho)+v1*np.sin(theta1-rho)*np.sin(rho+np.pi/2))
        
        v2x=((v2*np.cos(theta2-rho)*(m2-m1)+2*m1*v1*np.cos(theta1-rho))/
            (m2+m1)*np.cos(rho)+v2*np.sin(theta2-rho)*np.cos(rho+np.pi/2))
        v2y=((v2*np.cos(theta2-rho)*(m2-m1)+2*m1*v1*np.cos(theta1-rho))/
            (m2+m1)*np.sin(rho)+v2*np.sin(theta2-rho)*np.sin(rho+np.pi/2))
        
        body1.velocity.x=v1x
        body1.velocity.y=v1y
        
        body2.velocity.x=v2x
        body2.velocity.y=v2y
        
    def drainEnergy(self):
        ''' Drains energy (actually veclocity) from anything colliding with
        the detection box.
        '''
        # Set up collide algorithm
        circleCollide = pg.sprite.collide_circle_ratio(np.sqrt(2)/2.)
        
        # Build group of sprites on which to test collisions
        spriteGroup = pg.sprite.Group(self.listOfSprites[0:-1])
        col = pg.sprite.spritecollide(self.listOfSprites[-1],spriteGroup,
                                      dokill=False,collided=circleCollide)
        
        # If collision is true, drain energy
        for j in range(0,len(self.listOfSprites)):
            if self.listOfSprites[j] in col:
                
#=============================================================================
#                m = self.listOfSprites[j].body.mass
#                suckingParameter=1*10**23
#                self.listOfSprites[j].body.velocity.y *= (m-suckingParameter)/m
#                self.listOfSprites[j].body.velocity.x *= (m-suckingParameter)/m
#=============================================================================

                self.listOfSprites[j].body.velocity.y *= .99
                self.listOfSprites[j].body.velocity.x *= .99  
                
    def genSpriteAndBod(self, smallM, e, centralM, a, imageString, 
                        moon, keepCentral, omega=None, i=None):
        ''' Will add one, two, or three bodies and sprites to self.model.bodies 
        and self.listOfSprites. 
        
        Parameters
        -----
        smallM : list of floats
            A list consisting of the mass of the planet, and optionally the 
            mass of the moon.
        e : list of floats
            A list consisting of the desired eccentricities of the orbits.
        centralM : float
            The central mass of the system.
        a : list of floats
            The semimajor axis of the smaller bodies, first for the planet,
            second optionally for the moon.
        imageString : list of strings
            The strings of the file names of the images used as sprites
        moon : Boolean
            If moon==True, the function will look for the information for the 
            moon as the second etry of each list above.
            If moon==False, no moon will be added, and each list only needs one 
            entry.
        keepCentral : Boolean
            If keepCentral==True, the central reference body will be created
            and take the necessary information from the imageString list.
        omega : list of floats
            A list consisting of the desired angles for roation of the orbits.
            Defaults to zero.
        i : list of floats
            A list consisting of the desired angles of periasteron.
            Defaults to zero.
        '''
        # Check i and omega
        if i == None:
            i=np.zeros(len(smallM))
        if omega == None:
            omega=np.zeros(len(smallM))
            
        # Builds the base system
        bodies0=Simulation.SolarModel.closed2BodyAutoBuilder(omega[0], i[0], 
                                            centralM, smallM[0], a[0], e[0])
        bod00=bodies0[0]
        bod01=bodies0[1]
        
        if moon==True: # Builds a second system within the first
            
            bodies1=Simulation.SolarModel.closed2BodyAutoBuilder(omega[1], 
                                        i[1], smallM[0], smallM[1], a[1], e[1])
            bod10=bodies1[0]
            bod11=bodies1[1]

            # Reconcile orbits so that the first and second systems combine                        
            bod10.velocity=bod10.velocity+bod00.velocity
            bod10.position=bod10.position+bod00.position
            bod11.velocity=bod11.velocity+bod00.velocity
            bod11.position=bod11.position+bod00.position
            
            # Append bodies to the model's list of bodies
            self.model.bodies.append(bod10)
            self.model.bodies.append(bod11)
            
            # Load images and assemble sprites, append sprites
            bod10Image=pg.image.load(imageString[0]).convert_alpha()
            bod11Image=pg.image.load(imageString[1]).convert_alpha()
            
            self.listOfSprites.append(PlanetSprite(self.screen_size,
                                                   bod10Image,bod10,
                                                   self.scale))
            self.listOfSprites.append(PlanetSprite(self.screen_size,
                                                   bod11Image,bod11,
                                                   self.scale))
            
            if keepCentral == True : # Repeat steps to append the central body
                
                self.model.bodies.append(bod01)
                bod01Image=pg.image.load(imageString[2]).convert_alpha()
                self.listOfSprites.append(PlanetSprite(self.screen_size,
                                                       bod01Image,bod01,
                                                       self.scale))
            
        elif moon==False: # Skips building the second system
            
            # Load image, append body and sprites
            bod00Image=pg.image.load(imageString[0]).convert_alpha()
            self.model.bodies.append(bod10)
            self.listOfSprites.append(PlanetSprite(self.screen_size,
                                                   bod00Image,bod10,
                                                   self.scale))
            
            if keepCentral == True : # Same as above
                
                bod01Image=pg.image.load(imageString[1]).convert_apha()
                self.model.bodies.append(bod11)
                self.listOfSprites.append(PlanetSprite(self.screen_size,
                                                       bod11Image,bod11,
                                                       self.scale))
             
                
    def genSpriteAndBodManual(self, xPos, yPos, zPos, 
                              xVel, yVel, zVel, mass, imageString):
        ''' Generates new bodies and sprites with a manual set of parameters.
        
        Parameters
        -----
        xPos : list of float
            A list of all the x position values of the bodies being created.
        yPos : list of float
            A list of all the y position values of the bodies being created.
        zPos : list of float
            A list of all the z position values of the bodies being created.
        xPos : list of float
            A list of all the x velocity values of the bodies being created.
        yPos : list of float
            A list of all the y velocity values of the bodies being created.
        zPos : list of float
            A list of all the z velocity values of the bodies being created.
        mass : list of float
            A list of all the masses for each body being created
        imageString : list of strings
            A list of all the strings of the image file names for each body.
        '''
        
        bodies = Simulation.SolarModel.manualOrbitBuilder(xPos, yPos, zPos,
                                                     xVel, yVel, zVel, mass)
        for i in range(len(imageString)):
            bodImage=pg.image.load(imageString[i]).convert_alpha()
            self.listOfSprites.append(PlanetSprite(self.screen_size,bodImage,
                                                   bodies[i],self.scale))
            self.model.bodies.append(bodies[i])
            
    def initializeSimulation(self):
        ''' Runs the basic necessities to get the animation up and running.
        
        Parameters
        -----
        stepsize : float
            The stepsize used in the simulation. This gets overwritten, so this
            value is more just for initilizing than for actually determining
            how the simulation will run.
        secondsPerTick : float
            Actually determines how the system will run. Larger values means
            faster animation. Keep in mind this can be controlled during 
            simulation with the keypad + and - keys.
        '''
        # Miscellaneous bits
        stepsize=100000 # Doesn't do much
        secondsPerTick=self.sPT 
        screenInfo = pg.display.Info() # For making full screen work
        self.screen_size = (screenInfo.current_w,screenInfo.current_h) 
        self.t=[0]
        self.step = 0 
        self.detectionBox=False
        self.thrusting=False
        self.n=0
        self.BGcoords=[0,0]
        self.BGpos=[0,0]
        
        # Initialize simulation shtuff
        solver=Solver.RK4(stepsize,Physics.nBody.diffEq)
#        physics=Physics.UniformGravity(solver)
        physics=Physics.nBody(solver)
        stopCondition=StopConditions.afterTimeSetup(time=stepsize)     
        self.model=Simulation.SolarModel(stepsize,solver,physics,stopCondition)
        self.model.bodies=[]
        self.listOfSprites=[]
        
        # Initialize display surface
        self.screen = pg.display.set_mode([0,0], pg.FULLSCREEN)
        self.background = pg.image.load("starField2.jpg")
        self.screen.blit(self.background,(0,0))
        pg.display.update()
#        print pg.Surface.get_bytesize(self.background)
        
        # Initialize game over screen and win screen        
        self.GOscreen = pg.image.load("gameOver.jpg")
        self.GOscreen = pg.transform.scale(self.GOscreen, (self.screen_size[0],
                                                          self.screen_size[1]))
        self.gameOver = False
        self.winScreen = pg.image.load("trumpWinner.jpg")
        self.winScreen = pg.transform.scale(self.winScreen,
                                    (self.screen_size[0],self.screen_size[1]))
        self.win=False

        # Load other images
        self.redCrosshairs=pg.image.load('crosshairs100red.png')
        self.whiteCrosshairs=pg.image.load('crosshairs100.png')

        # Set up the timer.
        self.clock = pg.time.Clock()        
        self.ticks = pg.time.get_ticks()
        self.secondsPerTick = secondsPerTick
        
        # Load Sounds
#        pg.mixer.music.load("lowCarpet.wav")
#        pg.mixer.music.play(-1)
        self.drainSound=pg.mixer.Sound("bassHum.wav")
        self.thrustSound=pg.mixer.Sound("rocketThrust.wav")
        self.loseSound=pg.mixer.Sound("zenGong.wav")
#        self.frogs=pg.mixer.Sound("frogs.flac")
        
    def makeDetectionBox(self,image_file):
        ''' Makes a box that will be used in drainEnergy.
        
        Parameters
        -----
        image_file : pygame loaded image file
            The image for the detection box.
        '''
        box = PlanetSprite(self.screen_size,image_file,body=None,
                           scale=self.scale)
        self.listOfSprites.append(box)
        self.detectionBox=True # Really only used so as to not double up on
        # detecting collisions in both drainEnergy and the standard collision
        # detection method
        
    def orbVel(self,m1,m2,r):
        ''' Automatically calculates orbital velocities needed to make a 
        circular orbit around the central body. This velocity can easily be
        adjusted to get more elliptical orbits.
        
        Parameters
        -----
        m1 : float
            Mass of one body (usually the central mass).
        m2 : float
            Mass of the other body (usually the mass of the orbiting body).
        r : float
            The distance between the orbiting body and the center of mass.
        
        '''
        r=abs(r)
        G=6.6408*10**-11
        return (G*(m1+m2)/r)**0.5
        
    
    def processEvents(self):
        ''' Processes all user inputs.
        '''
        ship = self.listOfSprites[0] # This gets used a lot
        # Process events
        for event in pg.event.get():
            if event.type == pg.QUIT:
                self.done = True
            
            # Check for buttons that have been pressed
            elif event.type == pg.KEYDOWN:       
                # Escape quits
                if event.key == pg.K_ESCAPE:
                    self.done = True
                # '+' and '-' on the key pad speed up or slow down simulation
                if event.key == pg.K_KP_PLUS:
                    self.secondsPerTick *= 2
                if event.key == pg.K_KP_MINUS:
                    self.secondsPerTick *= 0.5 
                # '<' and '>' slow and speed up your thruster respectively
                if event.key == pg.K_PERIOD:
                    ship.thrust *=1.5
                if event.key == pg.K_COMMA:
                    ship.thrust *=.75
                # Up arrow gives thrust, plays sound, changes sprite
                if event.key == pg.K_UP:
                    ship.image=self.rot_center(self.rocketOnImageOrig,
                                               ship.angle)
                    self.thrusting=True
                    self.thrustSound.play()
                # Spacebar activates press plays sound and changes crosshairs
                if event.key == pg.K_SPACE:
                    self.drainSound.play()
                    self.listOfSprites[-1].image=self.redCrosshairs
                    
            # Check for buttons that have been released
            elif event.type == pg.KEYUP:
                # Releasing up arrow returns sprites and sounds to normal
                if event.key == pg.K_UP:
                    ship.image=self.rot_center(
                        self.rocketImageOrig,ship.angle)
                    self.thrusting=False
                    self.thrustSound.stop()
                # Releasing spacebar returns sound and crosshairs to normal
                if event.key == pg.K_SPACE:
                    self.drainSound.stop()
                    self.listOfSprites[-1].image=self.whiteCrosshairs

        # Process held buttons
        pressed = pg.key.get_pressed()
        # Left arrow rotates ship CCW
        if pressed[pg.K_LEFT]:
            if self.thrusting == True: # Determine which image to rotate
                img=self.rocketOnImageOrig
            else: img=self.rocketImageOrig
            ship.angle+=5 # Increase rotation by 5 degrees
            ship.angle%=360 # Check if divisible by 360
            # Rotate image about its center (doesn't actually work 100%)
            ship.image=self.rot_center(img,ship.angle)
        # Right arrow does all the same in the opposite direction
        if pressed[pg.K_RIGHT]:
            if self.thrusting == True:
                img=self.rocketOnImageOrig
            else: img=self.rocketImageOrig
            ship.angle-=5
            ship.angle%=360
            ship.image=self.rot_center(img,ship.angle)
        # Adds thrust in the appropriate amount on the appropriate axes
        if pressed[pg.K_UP]:
            yvel=np.cos(ship.angle*2*np.pi/360)*ship.thrust
            xvel=np.sin(ship.angle*2*np.pi/360)*ship.thrust
            self.model.bodies[0].velocity.y -=yvel
            self.model.bodies[0].velocity.x -=xvel
        # Down arrow gives a small amount of reverse thrust
        if pressed[pg.K_DOWN]:
            yvel=np.cos(ship.angle*2*np.pi/360)*ship.thrust/4.
            xvel=np.sin(ship.angle*2*np.pi/360)*ship.thrust/4.
            self.model.bodies[0].velocity.y +=yvel
            self.model.bodies[0].velocity.x +=xvel
        # Activate drain while spacebar is held
        if pressed[pg.K_SPACE]:
            self.drainEnergy()
            
    def rot_center(self, img, angle):
        '''Rotates a pygame surface, maintaining position. Doesn't actually 
        rotate about the center, unfortunately.
        
        Parameters
        -----
        img : pygame image
            The image to be rotated
        angle : float
            The number of degrees for the rotation CCW.
        '''
        loc = img.get_rect().center
        rot_sprite = pg.transform.rotate(img, angle)
        rot_sprite.get_rect().center = loc
        return rot_sprite
        
    def run(self):
        ''' The foundation for running the animation and simulation.
        '''
        self.done=False
        while not self.done:
            # Check user inputs
            self.processEvents()
            
            # Run collision detection
            if len(self.t)>1:
                self.checkCollisions()
#                print self.listOfSprites[0].body.position.y
#                self.checkFloorCollisions()
            
            # Check for win condition
            self.checkWin()
            
            # Check to see if the game has ended with a win or loss
            if self.gameOver != True and self.win != True: self.advance()
            
            # In case of a game over, play and blit appropriate things
            elif self.gameOver==True: 
                self.n+=1
                if self.n==1:
                    self.loseSound.play()
                    self.screen.blit(self.GOscreen,(0,0))
                    pg.display.update() 
                    
            # In case of a win, play and blit appropriate things
            elif self.win == True:
                self.n+=1
                if self.n==1:
#                    self.frogs.play()
                    self.screen.blit(self.winScreen,(0,0))
                    pg.display.update()
                    
        # If done, quit
        pg.quit()

        

#============================================================================#
#============================================================================#
#============================================================================#

class PlanetSprite(pg.sprite.Sprite):
    ''' Class for building and maintaining information about sprites.
    '''
    def __init__(self,screen_size,image_file,body,scale):
        
        pg.sprite.Sprite.__init__(self)
        self.body=body
        self.image = image_file
        self.rect = self.image.get_rect()
        self.screen_size=screen_size
        
        # Default assumptions
        self.blackHole=False
        self.mainBody =False
        self.angle=0
        self.thrust=0
        
        # Scale for drawing; the larger the absolute value of the exponent, the
        # closer together objects will appear.
        self.scale=scale
        
    def update(self, coords, pos, dim, angle):
        ''' Updates sprites based on their type and how you want things to
        move. Contains options for focusing on ship or on Center of Mass.
        '''
        
        self.coords=coords

        if self.body != None:           # Check to see if sprite has a body
            self.shipCenter(coords,pos,dim) # Use this to keep ship at center
#            self.centerOfMass()        # Use this to keep CoM at center
#            self.centerOfMass2()        # Use this to keep CoM at center

            
        else: # Used to manage detection box
#            print 'There is a detection box?'
            self.trackShip(angle)

    def shipCenter(self, coords, pos, dim):
        ''' Updates sprite based on the ship being at the center.
        '''
        body=self.body
        scale=self.scale
        if self.mainBody==False:
            self.rect.centerx = (body.position.x-pos[0])*scale + coords[0]+dim[0]/2.
            self.rect.centery = (body.position.y-pos[1])*scale + coords[1]+dim[1]/2.
        elif self.mainBody==True:
            self.rect.centerx = self.screen_size[0]/2.
#            print self.screen_size[0]/2.
            self.rect.centery = self.screen_size[1]/2.
#            print self.screen_size[1]/2.

    def centerOfMass(self):
        ''' Updates sprite based on the center of mass being at the center.
        '''
        body=self.body
        scale=self.scale
        self.rect.centerx = body.position.x*scale + self.screen_size[0]/2.
        self.rect.centery = self.screen_size[1] - body.position.y*scale - self.screen_size[1]/2.
        
    def centerOfMass2(self):
        ''' Updates sprite based on the center of mass being at the center.
        '''
        body=self.body
        scale=self.scale
#        print self.coords
        self.rect.centerx = (body.position.x*scale + self.screen_size[0]/2.)-self.coords[0]
        self.rect.centery = (self.screen_size[1]/2. - body.position.y*scale)-self.coords[1]
    
    def trackShip(self,angle):
        ''' Updates the detection box to follow the ship at the center.
        '''
        x=0
        y=self.image.get_rect()[3]/2.+20
        
        mag=y
        x=np.sin(angle/360.*2*np.pi)*mag
        y=np.cos(angle/360.*2*np.pi)*mag
        self.rect.centerx = self.screen_size[0]/2.-x
        self.rect.centery = self.screen_size[1]/2.-y
        
if __name__=='__main__':

    app = App(level)


    app.run()
        