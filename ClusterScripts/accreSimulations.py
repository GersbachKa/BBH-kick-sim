import numpy as np

import sys
sys.path.append('../')
from BBH_kick_sim import Simulator

import pickle

print('Done with imports')

if len(sys.argv)>1:
    args = sys.argv[1:]
else:
    args = []

#Command line arguments are in this order
#python accreSimulations.py {number of simulations} {collision radius} {use mass lookup} {random spin type:[uniform or zero]} {sortMass:[boolean]} {outputfile}

#e.g.
#python accreSimulations.py 10 2 False zero False /home/gersbaka/BBH-population-kicks/ClusterScripts/testout.pkl

params = {'cluster_mass': 1e6, 'radius': 2, 'imf_alpha':2.35, 'min_bh_star':10,
          'bh_mass_frac':0.5, 'use_mass_lookup':False, 'min_star':0.8, 'max_star':100,
          'vel_thresh':0.1,'rand_spin_type':'uniform'}

numSims                   =int(args[0])
params['radius']          =float(args[1])
params['use_mass_lookup'] = args[2] == 'True'
params['rand_spin_type']  =args[3]
sortMass                  = args[4] == 'True'
outfile                   =args[5]
print(sortMass)


seeds = []
all_gcs = []
for i in range(numSims):
    print(f'Running sim {i+1} of {numSims}')
    sim = Simulator.Simulator(params)
    seeds.append(sim.seed)
    sim.begin_sim(sort_mass_first=sortMass)
    all_gcs.append(sim.GC)


with open(outfile,'wb') as f:
    pickle.dump((params,seeds,all_gcs),f)
    
print('finished')