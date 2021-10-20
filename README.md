# Debris-covered glacier model

This repository contains files to run a full Stokes 2-D (flowline) ice-flow model of an example debris-covered glacier. A thorough description of the model and application results can be found in *Mayer C. and Licciulli C. (2021). The Concept of Steady State, Cyclicity and Debris Unloading of Debris-Covered Glaciers. Front. Earth Sci. 9:710276.* [doi: 10.3389/feart.2021.710276](https://doi.org/10.3389/feart.2021.710276).

The model is implemented using the Finite Element software [Elmer/Ice](http://elmerice.elmerfem.org/). Results in the paper were calculated using the Elmer version 8.3 (Rev: 05a2ed4). More recent Elmer versions should work as well (see Elmer [git repository](https://github.com/ElmerCSC/elmerfem)).

The most relevant file is the user function *USF_DebrisCoverage.f90*, which contains the implementation of debris advection, debris discharge mechanisms and the debris-influenced mass balance model.


## .sif files

* *Stokes_prognostic_no_debris.sif*: debris-free glacier evolution starting from an ice-free bedrock.
* *Stokes_prognostic.sif*: debris-induced glacier evolution starting from an arbitrary glacier state. Modify (and compile) the user function USF_DebrisCoverage.f90 in order to tune parameters governing debris dynamics and mass balance.

## Solvers

* *FreeSurfaceSolver_trans_stab.F90*: identical to *FreeSurfaceSolver.F90* from the Elmer distribution, but with addition of transient stabilization. Command to compile: ``elmerf90 FreeSurfaceSolver_trans_stab.F90 -o FreeSurfaceSolver_trans_stab.so``.

## USFs

* *USF_DebrisCoverage.f90*: user function to calculate debris advection, debris discharge and debris-influenced mass balance. Command to compile: ``elmerf90 USF_DebrisCoverage.f90 -o USF_DebrisCoverage``.

### Implemented discharge mechanisms:

* *DebrisSlideNode*: move entire debris column to the next sub-node, if surface slope exceeds a threshold.
* *DebrisSlideFull*: remove entire debris coverage downstream of a sub-node, if surface slope exceeds a threshold.
* *DebrisSlideYieldStress*: remove entire debris coverage downstream of a sub-node, if debris driving stress exceeds a threshold.
* *DebrisSlidePartial*: deterministic/random partial displacement of the debris column to the next sub-node combined with abrupt debris removal for high driving stresses.

Glacier systems without debris discharge can be reproduced setting all mechanisms to "*false*". We recommend to switch on discharge mechanisms only separately. The combination of more discharge mechanisms requires extensive testing and likely further coding.

### Relevant parameters:

* *DebrisEisRatio*: ratio between the height of the debris column within the debris/ice mixture and the total height of the debris/ice mixture.
* *MinHeigh*: minimum glacier thickness, below which nodes are considered ice-free.
* *MaxSlope*: maximum surface slope for the calculation of debris discharge.
* *YieldStress*: maximum debris driving stress for the calculation of debris discharge.
* *MonitoringFile*: file registering the total debris amount after specific calculation steps. In the current version the file must be deleted before restarting a new simulation.


