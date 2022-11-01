# BBH-population-kicks
The final project for Astr8020 for Kyle and Levi

Link to paper draft: https://www.overleaf.com/read/qdtkhqfzjdmx

# Research question

Given that we know black hole mergers can often lead to a velocity 'kick' after the merger, we want to know how important this effect is when it comes 
to a heirarchy of mergers. Specifically we want to know what happens to the population of black holes after subsequent mergers inside a globular 
cluster, a common place for black holes to form and grow.


## Research plan

To do this, we plan to interact populations of black holes with each other in subsequent mergers. Using a distribution of Masses, spins, and 
velocities, we will be able to semi-randomly interact black holes with each other using the SurfinBH python package, which will give us information 
about the resulting kick velocities, masses, and spins. Mass ratio information could also be recorded to verify our simulations are consistent with observational rates. The resulting population will be our second generation of black holes. Any black holes that leave the globular 
cluster will be recorded, but not interacted again. Each new generation will need to time evolve, using dynamical friction to slow velocities, we will 
re-interact these second generation black holes to make a third generation. Repeating this until we have some resulting population of intermediate mass 
black holes (i.e. 500 M_sun). This could also allow for inferences about the expected mass ratios for IMBH mergers.


Once we have many generations of black holes, we can bring the population densities together to get a total BBH population. 


To tie this back to LIGO, we can use the rates of BBH mergers as a way to change our percentage densities into physical densities.


### Some papers

- [Gravitational Wave Recoil and the Retention of Intermediate-Mass Black Holes](https://iopscience.iop.org/article/10.1086/591218)
This does almost everything we plan to do, though much more analytical/mathematical than what we intend to do.

- [SPIN-PRECESSION: BREAKING THE BLACK HOLE–NEUTRON STAR DEGENERACY](https://iopscience.iop.org/article/10.1088/2041-8205/798/1/L17)
- [Cosmological Black Hole Spin Evolution by Mergers and Accretion](https://iopscience.iop.org/article/10.1086/590379/meta)
These might be helpful when deriving initial spin parameter distributions

- [Global Stellar Budget for LIGO Black Holes](https://ui.adsabs.harvard.edu/abs/2020ApJ...889L..35J/abstract)
Helpful for black hole formation and merger rates

- [Large Merger Recoils and Spin Flips from Generic Black Hole Binaries](https://iopscience.iop.org/article/10.1086/516712/pdf)
Explores mass ratio 2 and highly misaligned spin mergers and resulting recoil

- [Binary black hole mergers from globular clusters: Masses, merger rates, and the impact of stellar evolution](https://journals.aps.org/prd/pdf/10.1103/PhysRevD.93.084029)
Similar study to ours but more extensive globular cluster modeling from 2016 using very few GW events.



