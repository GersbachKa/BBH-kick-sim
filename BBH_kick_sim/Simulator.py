import numpy as np
import sympy
import surfinBH
from copy import deepcopy
import time

from .GlobularCluster import GlobularCluster 
from .BlackHole import BlackHole


class Simulator:


    def __init__(self,params=None,print_missing=True,rand_seed=None):
        '''Params is a dictionary of values for the Globular Cluster Simulation.
        
        Params needs the following parameters, if any are left as none, defaults are provided
        
        cluster_mass - The total mass of the globular cluster in solar masses [1e6]
        radius - The radius where the BH collisions happen in parsecs [2]
        imf_alpha - The alpha value for the imf dn/dm = m**(-alpha) [2.35] 
        
        min_bh_star - Minimum mass of a star that forms a black hole in solar masses [10]
        bh_mass_frac - The fraction of the star mass that remains in the black hole [0.5]
        use_mass_lookup - Toggle whether to use the BH mass transform from Spera & Mapelli, 2017 [True] 
        
        min_star - The smallest star mass in the cluster in solar masses [0.8]
        max_star - The largest star mass in the cluster in solar masses [100]
         
        vel_thresh - Threshold velocity to allow a collision of black holes [0.1]
        rand_spin_type - The type of random for the random spin generator [uniform]
        '''
        if params==None:
            params = dict()
        pars = deepcopy(params)
        
        if 'cluster_mass' not in params.keys():
            pars.update({'cluster_mass':1e6})
            if print_missing:
                print(f"'cluster_mass' not set, defaulting to {pars['cluster_mass']}")
        if 'radius' not in params.keys():
            pars.update({'radius':2})
            if print_missing:
                print(f"'radius' not set, defaulting to {pars['radius']}")
        if 'imf_alpha' not in params.keys():
            pars.update({'imf_alpha':2.35})
            if print_missing:
                print(f"'imf_alpha' not set, defaulting to {pars['imf_alpha']}")
        if 'min_bh_star' not in params.keys():
            pars.update({'min_bh_star':10})
            if print_missing:
                print(f"'min_bh_star' not set, defaulting to {pars['min_bh_star']}")
        if 'bh_mass_frac' not in params.keys():
            pars.update({'bh_mass_frac':0.5})
            if print_missing:
                print(f"'bh_mass_frac' not set, defaulting to {pars['bh_mass_frac']}")
        if 'min_star' not in params.keys():
            pars.update({'min_star':0.8})
            if print_missing:
                print(f"'min_star' not set, defaulting to {pars['min_star']}")
        if 'max_star' not in params.keys():
            pars.update({'max_star':100})
            if print_missing:
                print(f"'max_star' not set, defaulting to {pars['max_star']}")
        if 'use_mass_lookup' not in params.keys():
            pars.update({'use_mass_lookup':False})
            if print_missing:
                print(f"'use_mass_lookup' not set, defaulting to {pars['use_mass_lookup']}")
        if 'vel_thresh' not in params.keys():
            pars.update({'vel_thresh':0.1})
            if print_missing:
                print(f"'vel_thresh' not set, defaulting to {pars['vel_thresh']}")
        if 'rand_spin_type' not in params.keys():
            pars.update({'rand_spin_type':'uniform'})
            if print_missing:
                print(f"'rand_spin_type' not set, defaulting to {pars['rand_spin_type']}")
                
        if pars['use_mass_lookup']:
            pars.update({'min_bh_star':16.2115})
            print('Using mass lookup table. Setting the minimum BH star to 16.2115')
            
        self.sim_params = pars
            
        if rand_seed==None:
            self.seed = np.random.randint(0,1_000_000_000)
        else:
            self.seed = int(rand_seed)
        
        self.rng = np.random.default_rng(self.seed)
        
        
        self.GC = GlobularCluster(pars['cluster_mass'],pars['radius'],pars['imf_alpha'],
                                  pars['min_bh_star'],pars['bh_mass_frac'],
                                  pars['min_star'],pars['max_star']
                                 )
        
        self._set_mass_transform(pars['imf_alpha'],pars['min_bh_star'],
                                 pars['bh_mass_frac'],pars['max_star'],pars['use_mass_lookup'])
        
        
        self.fit  = surfinBH.LoadFits('NRSur7dq4Remnant')
        self.fit2 = surfinBH.LoadFits('NRSur3dq8Remnant') #Used only with high q
        
        self.v_thresh = pars['vel_thresh']
        
        self.rand_spin_type = pars['rand_spin_type']
        
        for i in range(self.GC.Nbh):
            self.add_random_BH(id=str(i))
            
        self.collide_ratios = []
        print(f'Setup complete, Globular Cluster now has {len(self.GC.BHs)} black holes.')
         
    
    def begin_sim(self,stopTime=None,progress=1000,sort_mass_first=False):
        '''If stopTime==None, Run until one BH left'''
        cputime_s = time.time()
        
        if stopTime==None:
            stopTime=float(np.inf)
            print('No stop time specified, Running until 1 or 0 black holes remain')
        
        #Sort by mass first (As high mass are first to die)
        if sort_mass_first:
            self.sort_by_m_mass()
        
        col_count = 0
        init_BH = len(self.GC.BHs)
        while(len(self.GC.BHs)>1):
            if col_count%progress==0:
                print(f'{len(self.GC.BHs)}/{init_BH} remaining')
            
            #sort the list!
            self.sort_by_time()
            if self.GC.BHs[1].t < stopTime:
                #Take first two black holes
                bh1 = self.GC.remove_BH(0)
                bh2 = self.GC.remove_BH(0)
            
                #Collide them together
                mf,sf,vf,t = self.collideBH(bh1,bh2)
            
                #GC will check for ejected black holes,
                BHf, ejected = self.GC.add_BH(mf,sf,vf,t,bh1=bh1,bh2=bh2,ret=True,id=bh1.id+'+'+bh2.id)
            
                if not ejected:
                    #Time evolve the final BH to threshold velocity
                    self.time_evolve(BHf)
                col_count+=1
                
        print(f'Finished. Total Simulation time: {t}\nTotal CPU time: {time.time()-cputime_s}')

   
    def time_evolve(self,bh):
        #TODO: Include dynamical friction here
        #For now, just decrease velocity linearly
        a = -0.0001
        totV = np.sqrt(np.sum(np.square(bh.v))) #Get velocity magnitude
        
        if self.v_thresh<totV:
            t = (self.v_thresh-totV)/a
            bh.t+=t
            bh.v = self.v_thresh 
        else:
            pass
            #No need to decrease speed! Already at treshold!
       
        
    def collideBH(self,bh1,bh2):
        #Second BH always has the longer time, (sorted by time)
        t = bh2.t
        bh1.t = t #Collision happened at the same time!
        
        mtot = bh1.m+bh2.m
        
        if bh1.m>=bh2.m:
            q = bh1.m/bh2.m
            s1 = bh1.s_mag * self._random_uniform_sphere()
            bh1.s = s1
            s2 = bh2.s_mag * self._random_uniform_sphere()
            bh2.s = s2
        else:
            q = bh2.m/bh1.m
            s1 = bh2.s_mag * self._random_uniform_sphere()
            bh2.s=s1
            s2 = bh1.s_mag * self._random_uniform_sphere()
            bh1.s=s2
        
        if q<6:
            mf,chif,vf,_,_,_ = self.fit.all(q,s1,s2,allow_extrap=True)
        else:
            #big mass ratio, use second fitter
            #print(f'Large mass ratio encountered: q={q}')
            s1 = np.sqrt(np.sum(np.square(s1)))*self.rng.choice([[0,0,1],[0,0,-1]])
            s2 = np.sqrt(np.sum(np.square(s2)))*self.rng.choice([[0,0,1],[0,0,-1]])
            mf,chif,vf,_,_,_ = self.fit2.all(q,s1,s2,allow_extrap=True)
        mf *= mtot
        vf *= 3e5
        self.collide_ratios.append(q)
        return mf,chif,vf,t
    
    
    def sort_by_time(self):
        '''Sort the Black holes by the time'''
        self.GC.BHs.sort(key=lambda b: b.t)
        
    def sort_by_m_mass(self):
        self.GC.BHs.sort(key=lambda b: -b.m)
    
    
    def add_random_BH(self,id=None):
        '''Make a random BH to add to the GC. Starts at time 0 with 0 velocity'''
        m = self.random_mass()
        s = self.random_spin()
        v = 0
        t = 0
        
        return self.GC.add_BH(m,s,v,t,bh1=None,bh2=None,id=id)
        
    
    def random_mass(self):
        return self._mass_transform(self.rng.random())
      
        
    def random_spin(self):
        '''SurfinBH works best with spin [0,0.8) but works 'well' up to 1'''
        if self.rand_spin_type=='uniform':
            return self.rng.random()*(0.9-0) * self._random_uniform_sphere()
        elif self.rand_spin_type=='zero':
            return 0 * self._random_uniform_sphere()
    
    
    def _set_mass_transform(self,imf,min_bh_star,bh_mass_frac,max_star,use_mass_lookup):
        '''sets up the random mass generator. Uses inverse transform sampling'''
        print('Setting up analytic mass distribution. '+ 
              'This may take a while depending on your imf alpha')
        
        m = sympy.Symbol('m',real=True)
        p = sympy.Symbol('p',real=True)
        func = m**(-(imf-1))
        
        min_bh_mass = min_bh_star
        max_bh_mass = max_star

        norm = sympy.integrate(func, (m, min_bh_mass, max_bh_mass))
        
        pdf = (1/norm)*m**(-(imf-1))
        cdf = sympy.integrate(pdf,(m,min_bh_star,m))
        
        quantile = sympy.solvers.solve(cdf-p,m)
        
        try:
            quantile = sympy.lambdify(p,quantile[0])
        except:
            quantile = sympy.lambdify(p,quantile)
                      
        if use_mass_lookup:
            func = self._make_mass_lookup()
            mt = lambda x: float(func(quantile(x)))
            self._mass_transform = mt
        else:
            mt = lambda x: bh_mass_frac*quantile(x)
            self._mass_transform = mt
        print('Done')
        
    
    def _random_uniform_sphere(self):
        '''Returns a random vector of a unit sphere in cartesian'''
        phi = self.rng.random()*2*np.pi
        theta = np.arcsin(self.rng.random())
        xyz = np.array([np.cos(phi)*np.sin(theta), np.sin(phi)*np.sin(theta), np.cos(theta)])
        return xyz
    
    def _make_mass_lookup(self):
        from scipy.interpolate import interp1d
                      
        massLookup = np.array([
            [8.0,1.3], 
            [13.8,1.4], 
            [19.5,4.0],
            [25.2,10.7], 
            [31.0,21.3], 
            [36.8,33.0], 
            [42.5,38.1], 
            [48.2,43.2], 
            [54.0,48.3],
            [59.8,53.4],
            [65.5,48.5],
            [71.2,28.8],
            [77.0,31.1],
            [82.8,33.8],
            [88.5,36.6],
            [94.2,39.2],
            [100.0,41.9]
        ])
        
        return interp1d(massLookup[:,0],massLookup[:,1])
        
        