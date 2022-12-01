
import numpy as np

import sys
sys.path.append('..')
from BBH_kick_sim import Simulator

import pickle

if len(sys.argv)>1:
    args = sys.argv[1:]
else:
    args = []

#Command line arguments are in this order
#python accreSimulations.py {number of simulations} {cluster mass} {collision radius} {imf} {minimum bh star mass} {bh mass fraction} {use mass lookup} {minimum star mass} {maximum star mass} {merger velocity threshold} {random spin type:[uniform or zero]} {sortMass:[boolean]} {outputfile}

#e.g.
#python accreSimulations.py 10 1e6 2 2.25 10 0.5 False 0.8 100 0.1 zero False /Users/gersbaka/Documents/Homework/Astro8020/BBH-population-kicks/ClusterScripts/testout.pkl

params = {'cluster_mass': 1e6, 'radius': 2, 'imf_alpha':3.25, 'min_bh_star':10,
          'bh_mass_frac':0.5, 'use_mass_lookup':False, 'min_star':0.8, 'max_star':100,
          'vel_thresh':0.1,'rand_spin_type':'uniform'}

numSims                   =int(args[0])
params['cluster_mass']    =float(args[1])
params['radius']          =float(args[2])
params['imf_alpha']       =float(args[3])
params['min_bh_star']     =float(args[4])
params['bh_mass_frac']    =float(args[5])
params['use_mass_lookup'] =bool(args[6])
params['min_star']        =float(args[7])
params['max_star']        =float(args[8])
params['vel_thresh']      =float(args[9])
params['rand_spin_type']  =args[10]
sortMass                  =bool(args[11])
outfile                   =args[12]


seeds = []
all_gcs = []
for i in range(numSims):
    sim = Simulator.Simulator(params)
    seeds.append(sim.seed)
    sim.begin_sim(sort_mass_first=sortMass,progress=10000)
    all_gcs.append(sim.GC)


with open(outfile,'wb') as f:
    pickle.dump((params,seeds,all_gcs),f)
    
print('finished')
    
