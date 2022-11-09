import numpy as np
import sympy
import surfinBH

class BlackHole:
    '''The Black hole object, used to keep track of their parameters'''
    def __init__(self,mass,spin,velocity,time,bh1=None,bh2=None):
        '''Creating a new Black hole object'''
        
        # Need to keep track of formation parameters and time evolved
        self.m = mass
        
        self.s = spin
        #self.s_i = spin #Does spin evolve over time?
        
        self.v = velocity
        self.v_i = velocity
        
        self.t = time
        self.t_i = time
        
        self.parents = (bh1,bh2)
        
    def add_velocity_time(velocity,time):
        self.v += velocity
        self.t += time
