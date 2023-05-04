#!/usr/bin/env python

from distutils.core import setup

setup(name='SpaceVehiclePropagator',
      version='1.0',
      description='A Kepler STM Space Vehicle Propagator',
      author='Dean Keithly',
      author_email='drk94@cornell.edu',
      url='https://github.com/deanthedream/SpaceVehiclePropagator',
      packages=['EarthOrbiter', 'Plotter', 'KOE'],
     )