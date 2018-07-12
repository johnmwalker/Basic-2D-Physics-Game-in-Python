# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 14:51:41 2016

@author: John
"""
import numpy as np

class ThermalBody(object):
    '''Creates a body object whose only attribute is a temperature.
    
    Attributes
    -----
    temperature : float
        The temperature of the body.
    '''
    def __init__(self, temperature):
        self.temperature = temperature
        
class GravBody(object):
    '''Creates a body object whose attributes are velocity and height.
    
    Attributes
    -----
    velocity : instance of Vector
        A vector containing the three velocity vectors, Vx, Vy, and Vz.
    position : instance of Vector
        A vector containing the three position vectors, x, y, and z (or height).
    mass : float
        Mass of body in kg.
    '''
    def __init__(self,velocity,position,mass,radius=None):

        self.position = position
        self.velocity = velocity
        self.mass = mass
        self.radius = radius
        
    def serialize(self):
        
        vx = self.velocity.x
        vy = self.velocity.y
        vz = self.velocity.z
        x = self.position.x
        y = self.position.y
        z = self.position.z
        
        return np.array([vx,vy,vz,x,y,z])
        
    def deserialize(self,array):
        
        self.velocity.x=array[0]
        self.velocity.y=array[1]
        self.velocity.z=array[2]
        self.position.x=array[3]
        self.position.y=array[4]
        self.position.z=array[5]
        