import numpy as np
import sympy
import surfinBH

from .BlackHole import BlackHole


class GlobularCluster:
    
    def __init__(self,mass,radius,
                 imf_alpha=2.35,min_bh_star=10,bh_mass_frac=0.5,
                 min_star=0.8,max_star=100):
        ''''''

        #IMF and BH information (For N_BH calculation)
        self.imf_alpha = imf_alpha
        self.m_bh_star_min = min_bh_star #In solar masses
        self.bh_mass_fraction = bh_mass_frac 
        
        #Min and max star masses (For N_BH calculation)
        self.m_star_min = min_star
        self.m_star_max = max_star
             
        #Set total mass and radius of the globular cluster
        self.mass = mass #GC_mass must be in solar masses
        self.radius = radius * 3.086e16 #Radius in PC. Convert to m
        
        #Calculate escape velocity
        G = 6.6743e-11
        m_sun = 2e30
        self.v_esc = np.sqrt((2*G*self.mass*m_sun)/self.radius) *(1e-3) #Extra factor for m/s -> km/s
        
        self.Nbh = self._calculateNbh()
        self.BHs = []
        self.ejected = []
        
    def _calculateNbh(self):
        '''Calculate number of BHs given a GC mass
        *Uses Salpeter IMF'''
        
        #Need to find normalizing constant
        m = sympy.Symbol('m',real=True)
        func = m**(-(self.imf_alpha-1))
        a = self.mass/sympy.integrate(func,(m,self.m_star_min,self.m_star_max))
        
        #Now use in number integral
        func2 = a*m**(-self.imf_alpha) 
        Nbh = sympy.integrate(func2,(m,self.m_bh_star_min,self.m_star_max))
        return int(Nbh)
    
    
    def add_BH(self,mass,spin,velocity,time,bh1=None,bh2=None,ret=False):
        bh = BlackHole(mass,spin,velocity,time,bh1,bh2)
        if np.sqrt(np.sum(np.square(velocity))) < self.v_esc:
            #Still in system
            eject = False
            self.BHs.append(bh)
            if ret:
                return self.BHs[-1], eject
        else:
            eject = True
            self.ejected.append(bh)
            if ret:
                return self.ejected[-1], eject
        
        
    def remove_BH(self,index):
        return self.BHs.pop(index)
        
