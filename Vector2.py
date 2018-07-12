# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 13:56:08 2016

@author: John
"""
import math
import numpy as np

class Vector(object):
    '''Generic Vector class
    
    Attributes
    -----
    x : float
        cartesian x component
        
    y : float
        cartesian y component
            
    r : float
        polar radial component
        
    theta : float
        polar azimuthal angle
    '''
    
    def __init__(self, x=None,y=None):
        ''' Let's get set up!'''
        
        self.r = math.sqrt(x**2+y**2)
        self.theta = math.atan2(y, x)
      
    @property
    def x(self):
        return math.cos(self.theta)*self.r
        
    @x.setter
    def x(self, x):

        r = self.r     
        y = self.y
        theta = self.theta
        
        self.r = math.sqrt(x**2 + y**2)
        self.theta = math.atan2(y, x)
        
    @property
    def y(self):
        return math.sin(self.theta)*self.r
        
    @y.setter
    def y(self, y):
        r = self.r
        x = self.x       
        theta = self.theta
        
        self.r = math.sqrt(x**2 + y**2)
        self.theta = math.atan2(y, x)
        
    def __repr__(self):
        return "Vector(r=%f, y=%f)"%(self.r,self.y)
        
def main():
    
    X = 1
    Y = 2    
    
    V = Vector(x=X, y=Y)
    
    if V.r != math.sqrt(X**2+Y**2):
        print 'You failed!'
    else: print 'Your radius is correct at ' + repr(V.r)
        
    if V.theta != math.atan2(Y,X):
        print 'You failed!'
    else: print 'Your theta is correct at ' + repr(V.theta)
    
    print "Now we're going to reassign r and theta!"
    
    R = 1
    Theta = math.pi    
    
    V.r = R
    V.theta = Theta
    
    if V.x != math.cos(Theta)*R:
        print 'You failed!'
    else: print 'Your x value is correct at ' + repr(V.x)
    
    if V.y != math.sin(Theta)*R:
        print 'You failed!'
    else: print 'Your y value is correct at ' + repr(V.y)
    
    return V
    
if __name__ == '__main__': V = main()