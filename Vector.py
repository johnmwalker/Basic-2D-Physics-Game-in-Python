# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 13:56:08 2016

@author: John
"""
import math
import numpy as np

class Vector(object):
    '''Generic Vector class in either two or three dimensions.
    
    So far supports basic boolean == or !=, addition, multiplication, and
    subtraction between Vector objects. 
    
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
    
    def __init__(self, x=None,y=None,z=None,r=None,theta=None,phi=None):
        ''' Let's get set up!'''
        if x!=None:
            self.x = x
            self.y = y
            self.z = z
        else:
            self.x = r*math.sin(phi)*math.cos(theta)
            self.y = r*math.sin(phi)*math.sin(theta)
            self.z = r*math.cos(phi)
            
    def rotX(self,theta):
        '''
        '''
        row1=np.array([1,0,0])
        row2=np.array([0,math.cos(theta),math.sin(theta)])
        row3=np.array([0,-math.sin(theta),math.cos(theta)])
        rotMatrix=np.vstack([row1,row2,row3])
        
        xyzMatrix=np.array([self.x,self.y,self.z])
        newxyzMatrix=np.dot(xyzMatrix,rotMatrix)
#        newxyzMatrix=np.sum(newxyzMatrix,axis=0)
        self.x=newxyzMatrix[0]
        self.y=newxyzMatrix[1]
        self.z=newxyzMatrix[2]
    
    def rotY(self,theta):
        '''
        '''
        row1=np.array([math.cos(theta),0,-math.sin(theta)])
        row2=np.array([0,1,0])
        row3=np.array([math.sin(theta),0,math.cos(theta)])
        rotMatrix=np.vstack([row1,row2,row3])
        
        xyzMatrix=np.array([self.x,self.y,self.z])
        newxyzMatrix=np.dot(xyzMatrix,rotMatrix)
#        newxyzMatrix=np.sum(newxyzMatrix,axis=0)
        self.x=newxyzMatrix[0]
        self.y=newxyzMatrix[1]
        self.z=newxyzMatrix[2]
        
    def rotZ(self,omega):
        '''
        '''
        row1=np.array([math.cos(omega),math.sin(omega),0])
        row2=np.array([-math.sin(omega),math.cos(omega),0])
        row3=np.array([0,0,1])
        rotMatrix=np.vstack([row1,row2,row3])
        
        xyzMatrix=np.array([self.x,self.y,self.z])
        newxyzMatrix=np.dot(xyzMatrix,rotMatrix)
#        newxyzMatrix=np.sum(newxyzMatrix,axis=0)
        self.x=newxyzMatrix[0]
        self.y=newxyzMatrix[1]
        self.z=newxyzMatrix[2]
      
    @property
    def r(self):
        return math.sqrt(self.x**2+self.y**2+self.z**2)
        
    @r.setter
    def r(self, r):
        x = self.x
        y = self.y   
        z = self.z
        theta = self.theta
        phi = self.phi
        
        self.x = r*math.sin(phi)*math.cos(theta)
        self.y = r*math.sin(phi)*math.sin(theta)
        self.z = r*math.cos(phi)
        
    @property
    def theta(self):
        theta=math.atan2(self.y, self.x)
        if theta < 0:
            theta = 2*math.pi+theta
        return theta
        
    @theta.setter
    def theta(self, theta):
        x = self.x
        y = self.y        
        r = self.r
        phi = self.phi
        
        self.x = r*math.sin(phi)*math.cos(theta)
        self.y = r*math.sin(phi)*math.sin(theta)
        
    @property
    def phi(self):
        return math.atan2(math.sqrt(self.x**2+self.y**2),self.z)
        
    @phi.setter
    def phi(self, phi):
        x = self.x
        y = self.y  
        z = self.z
        r = self.r
        theta = self.theta
        
        self.x = r*math.sin(phi)*math.cos(theta)
        self.y = r*math.sin(phi)*math.sin(theta)
        self.z = r*math.cos(phi)
        self.r = math.sqrt(x**2+y**2+z**2)
        
#    @property
#    def x(self):
#        return self.r*math.sin(self.phi)*math.cos(self.theta)
#        
#    @x.setter
#    def x(self, x):
#
#        r = self.r     
#        y = self.y
#        z = self.z
#        theta = self.theta
#        
#        self.r = math.sqrt(x**2 + y**2 + z**2)
#        self.theta = math.atan2(y, x)
#        
#    @property
#    def y(self):
#        return self.r*math.sin(self.phi)*math.sin(self.theta)
#        
#    @y.setter
#    def y(self, y):
#        r = self.r
#        x = self.x       
#        theta = self.theta
#        
#        self.r = math.sqrt(x**2 + y**2 + z**2)
#        self.theta = math.atan2(y, x)
#        
#    @property
#    def z(self):
#        return self.r*math.cos(self.phi)
#        
#    @z.setter
#    def z(self, z):
#        r = self.r
#        x = self.x     
#        y = self.y
#        theta = self.theta
#        
#        self.r = math.sqrt(x**2 + y**2 + z**2)
#        self.phi = math.atan2(math.sqrt(x**2+y**2),z)

        
    def __repr__(self):
        return "Vector(x=%f, y=%f, z=%f)"%(self.x,self.y,self.z)
        
    def __mul__(self,v):
        return self.x*v.x+self.y*v.y+self.z*v.z

    def __add__(self,v):
        xNew = self.x+v.x
        yNew = self.y+v.y
        zNew = self.z+v.z
        return Vector(x=xNew, y=yNew, z=zNew)

    def __sub__(self,v):
        xNew = self.x-v.x
        yNew = self.y-v.y
        zNew = self.z-v.z
        return Vector(x=xNew, y=yNew, z=zNew)
    
#    def __lt__(self, other):
#        return self.pages < other
#    
#    def ___le__(self, other):
#        return self.pages <= other
    
    def __eq__(self, other):
        if self.x==other.x and self.y==other.y and self.z==other.z: 
            return True
        else:
            return False
    
    def __ne__(self, other):
        if self.x!=other.x or self.y!=other.y or self.z!=other.z:
            return True
        else:
            return False
    
#    def __gt__(self, other):
#        return self.pages > other
#    
#    def __ge__(self, other):
#        return self.pages >= other
    
    def __abs__(self):
        return Vector(x=abs(self.x),y=abs(self.y),z=abs(self.z))
        
#def main():
    
#    print '\n', "First tests are in two dimensions.", '\n'
#        
#    X = 1
#    Y = 2  
#    
#    V = Vector(x=X, y=Y)
#
####### Two-Dimensional Tests   
# 
#    if V.r != math.sqrt(X**2+Y**2):
#        print 'You failed to calculate radius!'
#    else: print 'Your radius is correct at ' + repr(V.r)
#        
#    if V.theta != math.atan2(Y,X):
#        print 'You failed to calculate theta!'
#    else: print 'Your theta is correct at ' + repr(V.theta)
#    
#    print "Now we're going to reassign r and theta!"
#    
#    R = 1
#    Theta = math.pi    
#    
#    V.r = R
#    V.theta = Theta
#    
#    if V.x != math.cos(Theta)*R:
#        print 'You failed to get x based on a new radius and theta!'
#    else: print 'Your x value is correct at ' + repr(V.x)
#    
#    if V.y != math.sin(Theta)*R:
#        print 'You failed to get a new y based on a new radius and theta!'
#    else: print 'Your y value is correct at ' + repr(V.y)
    
###### Three-Dimensional Tests

#    print('\n', "Next tests are in three dimensions.", '\n')
#
#    X = 1
#    Y = 2  
#    Z = 3
#    
#    V = Vector(x=X, y=Y, z=Z)
#    
#    if V.r != math.sqrt(X**2+Y**2+Z**2):
#        print('You failed to set radius based on x, y, and z!!!!!!')
#    else: print('Your radius is correct at ' + repr(V.r))
#        
#    if V.theta != math.atan2(Y,X):
#        print('You failed to set theta based on x, y, and z!!!!!!!')
#    else: print('Your theta is correct at ' + repr(V.theta))
#    
#    if V.phi != math.atan2(math.sqrt(X**2+Y**2),Z):
#        print('You failed to set phi based on x, y, and z!!!!!!!')
#    else: print('Your phi is correct at ' + repr(V.theta))
#    
#    print('\n', "Now we're going to reassign r, theta, and phi!")
#    
#    R = 10
#    Theta = math.pi   
#    Phi = math.pi/2.
#    
#    V.r = R
#    V.theta = Theta
#    V.phi = Phi
#    
#    if V.x != math.sin(Phi)*math.cos(Theta)*R:
#        print 'You failed to set a new x value based on a new radius and theta!!!!!'
#    else: print 'Your x value is correct at ' + repr(V.x)
#    
#    if V.y != math.sin(Phi)*math.sin(Theta)*R:
#        print 'You failed to set a new y value based on a new radius and theta!!!!!!'
#    else: print 'Your y value is correct at ' + repr(V.y)
#    
#    if V.z != math.cos(Phi)*R:
#        print 'You failed to set a new z value based on a new radius and theta!!!!!!'
#    else: print 'Your z value is correct at ' + repr(V.y)
#       
#    print '\n', "Now we're going to try some basic mathematical functions!"
#    
#    x1,y1,z1=(1,1,1)    
#    x2,y2,z2=(2,2,2)
#    
#    V1 = Vector(x=x1,y=y1,z=z1)
#    V2 = Vector(x=x2,y=y2,z=z2)
#    
#    if V1*V2 != x1*x2+y1*y2+z1*z2:
#        print 'You failed to multiply successfully!'
#    else: print 'Multiplication succeeded with the result ' + repr(V1*V2)
#    
#    if V1+V2 != Vector(x=x1+x2,y=y1+y2,z=z1+z2):
#        print 'You failed to add successfully!'
#    else: print 'Addition succeeded with the result ' + repr(V1+V2)
#    
#    if V1-V2 != Vector(x=x1-x2,y=y1-y2,z=z1-z2):
#        print 'You failed to subtract successfully!'
#    else: print 'Subtraction succeeded with the result ' + repr(V1-V2)
#
#    return V, V1, V2
    
#if __name__ == '__main__': V, V1, V2 = main()