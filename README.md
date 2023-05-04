# SpaceVehiclePropagator
A Two-Body Space Vehicle Propagator

First, this kind of code has been written hundreds to thousands of times and the most appropriate solution would be to use one that has already been coded and optimized.
My advisor from Cornell had https://github.com/dsavransky/keplertools which I used for propagating exoplanetary systems.
It has cython optimization and should be close to the fastest that can be achieved within python.

Since that didn't seem like the spirit of the request; I use the formulation from Algorithm 8 of Vallado for the propagation which uses the kepler state transition matrix for propagation.
I additionally included some plots that make use of basemap; which enables the plotting of an Earth.
Plotting like this typically works better in 3D.
I am not pleased that I was not able to get the axis ticks to display in the plot containing the Earth, but I am confident I could find a solution given more time.

In general, I only propagated a single planet and used/save intermediate perameters to avoid duplicate computation at the marginal cost of additional memeory.
I put the individual KOE into an EarthOrbiter object instead of using the numpy vector math to handle computations for multiple KOE.
There could be marginally better comments in a few of the functions. Lint rules were loosely followed.


If given more time, it would have been interesting to propagate for an orbit using true anomaly, create a fake TLE from the KOE, and formulate an optimization function with scipy wrapping an SPG4 propagator varying the TLE parameters.
This would find the error minimizing TLE and the SPG4 propagator can be used for propagation.



## Installation

Install git and use the github link to clone to a local repository
From the SpaceVehiclePropagator folder, run 
```
python setup.py install
```

Or with PYPI by
```
pip install SpaceVehiclePropagator
```

## Running
As an example, from the SpaceVehiclePropagator folder, run 
```
python ./EarthOrbiter/EarthOrbiter.py --a 6771 --e 0.02 --i 0 --W 0 --w 0 --v0 0 --dt 60 --mu 398600.4415
```

Alternatively, from an ipython session in the SpaceVehiclePropagator folder, you can run
```
from EarthOrbiter import EarthOrbiter
EarthOrbiter.main()
```


## Documentation

https://SpaceVehiclePropagator.readthedocs.io/
