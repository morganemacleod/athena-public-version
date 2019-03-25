athena
======
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

Athena++ radiation MHD code


this branch:
============

This branch includes algorithms to create modified, "merged", zones in the PHI direction near coordinate poles (theta=0, pi). Polar averaging averages conserved quanitites within sets of zones near the poles to create artifically enlarged zones with less extreme aspect ratios. The changes are within the hydro class and are mostly contained within the file "average_conserved.cpp". The timestep is relaxed accordingly.

In practice, in the zone *nearest* the pole, all x3-direction zones within a meshblock are averaged. From there averaging drops off with distance from the pole to preserve zone aspect ratios.  For example, if there are 128 zones in phi, and the meshblocks are 16x16x16, there will be 8 zones around the pole. 

There are two versions of how to do the momentum averaging in average_conserved.cpp, the first is a scalar average, the second is a vectorial average. We find that the vectorial average approach allows waves to pass more easily through the poles. It loses explicit conservation of momentum in the x3direction, however, so there may be situations where the scalar-average is preferable. 

A test problem called "pwave" is provided that sends a low-amplitude pressure wave through the pole. An input file is provided in the inputs/hydro directory. Users can toggle the "polar_average" parameter in the "hydro" block. 

Development of these techniques by Morgan MacLeod and Wenrui Xu. 
if you use in publications, please cite:
https://ui.adsabs.harvard.edu/#abs/2018ApJ...863....5M/abstract, https://ui.adsabs.harvard.edu/#abs/2018ApJ...868..136M/abstract, 
and W. Xu et al. (in prep, 2019) 
and please link to this repository. 

Changes, improvements, suggestions are welcomed! 
