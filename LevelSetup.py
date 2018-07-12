# -*- coding: utf-8 -*-
"""
Created on Tue May 17 19:39:11 2016

@author: John
"""
import Body
import Physics
import Search
import Simulation
import Solver
import Vector
import numpy as np
import math
import matplotlib.pyplot as plt
import StopConditions
import Error
import copy as cp
import pygame as pg
import TestStation as ts

class LevelSetup():
    '''
    '''
    def __init__(self, levelNumber):
        app=ts.App()
        if levelNumber==1: 
            return self.level1()
            
    def level1(self):
        ##### Parameters for Autobuilder #####

        M         =2.*10**34
        smallM    =[6.*10**30, 7.*10**22]
        a         =[4.0*10**12, 3.0*10**11]
        e         =[0,0]
        imageString = ["gasPlanetSmall.png", "ganymedeSmall.png"]
        
        ##### Parameters for Manual builder #####

        mass=[M, 1*10**23  , 1*10**24 , 1*10**25 , 1*10**26 ]
        xPos=[0, 1*10**12, -1*10**12, 2.*10**12, 3.*10**12]
        yPos=[0, 0,          0,         0,         0        ]
        zPos=[0, 0,          0,         0,         0        ]
        
        vel=app.orbVel(M,np.array(mass[1:]),np.array(xPos[1:]))
        
        xVel=[0, 0,          0,         0,         0        ]
        yVel=[0, vel[0],     -vel[1],    vel[2],    vel[3]   ]
        zVel=[0, 0,          0,         0,         0        ]
        imageString2 =["blackhole.png","planetSmall.png","planet2Small.png",
                                       "planet3Small.png","planet4Small.png"]
        
        ##### Build Rocket #####
        
        mRocket = [6.*10**30]
        xRocket,yRocket,zRocket    = [5.0*10**12], [0], [0]
        vxRocket,vyRocket,vzRocket = [0],[app.orbVel(M,mRocket[0],xRocket[0])],[0]
        rocketImage = ["shipSmall.png"]
        
        app.genSpriteAndBodManual(xRocket,yRocket,zRocket,
                                   vxRocket,vyRocket,vzRocket,
                                   mRocket,rocketImage)
                                   

                                   
        # Initialize ship characteristics
                                   
        self.listOfSprites[0].angle= 0
        self.listOfSprites[0].thrust = 20000.
        
        ##### Build bodies #####
        
        ts.App.genSpriteAndBodManual(xPos, yPos, zPos, xVel, yVel, zVel, mass, 
                                                               imageString2)  
        ts.App.genSpriteAndBod(smallM,e,M,a,imageString,
                        moon=True,keepCentral=False)


                        
        self.listOfSprites[1].blackHole = True
        self.listOfSprites[0].mainBody = True
        
        return self.listOfSprites
