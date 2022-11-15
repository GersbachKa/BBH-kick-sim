import numpy as np
import sympy
import surfinBH

class BlackHole:
    '''The Black hole object, used to keep track of their parameters'''
    def __init__(self,mass,spin,velocity,time,bh1=None,bh2=None):
        '''Creating a new Black hole object'''
        
        # Need to keep track of formation parameters and time evolved
        self.m = mass #Solar masses
        
        self.s = spin #Unitless spin
        self.s_mag = np.sum(np.square(spin))
        #self.s_i = spin #Does spin evolve over time?
        
        self.v = velocity #Km/s
        self.v_i = velocity
        
        self.t = time #s
        self.t_i = time
        
        self.parents = (bh1,bh2)
    
    def n_parents(self):
        n = 0
        if self.parents[0]==None:
            return n
        n += 1 + self.parents[0].n_parents()
        n += 1 + self.parents[1].n_parents()
        
        return n
    
    def __str__(self):
        BH = f'[M:{round(self.m,3)},S_mag:{round(self.s_mag,3)},'
        BH+= f'v_mag:{round(np.sqrt(np.sum(np.square(self.v))),3)},'
        BH+= f'v_imag:{round(np.sqrt(np.sum(np.square(self.v_i))),3)},'
        BH+= f't:{round(self.t,3)},t_i:{round(self.t_i,3)},'
        BH+= f'n_parents:{self.n_parents()}]'
        return BH
