# SpaceVehiclePropagator
A Two-Body Space Vehicle Propagator

First, this kind of code has been written hundreds to thousands of times and the most appropriate solution would be to use one that has already been coded and optimized.
My advisor from Cornell had https://github.com/dsavransky/keplertools which I used for propagating exoplanetary systems.
It has cython optimization and should be close to the fastest that can be achieved within python.

Since that didn't seem like the spirit of the request; I use the formulation from Algorithm 8 of Vallado for the propagation.
I additionally included some plots that make use of basemap; which enables the plotting of an Earth.
Plotting like this typically works better in 3D.
I am not pleased that I was not able to get the axis ticks to display in the plot containing the Earth, but I am confident I could find a solution given more time.





## Installation

Install git and use the github link to clone to a local repository
From the SpaceVehiclePropagator folder, run 
```
python setup.py install
```

## Running
From the SpaceVehiclePropagator folder, run 
```
python ./EarthOrbiter/EarthOrbiter.py
```


## Documentation

https://SpaceVehiclePropagator.readthedocs.io/
