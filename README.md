# SpaceVehiclePropagator
A Two-Body Space Vehicle Propagator





![Build Status](https://github.com/dsavransky/keplertools/actions/workflows/ci.yml/badge.svg)
[![Coverage Status](https://coveralls.io/repos/github/dsavransky/keplertools/badge.svg?branch=main)](https://coveralls.io/github/dsavransky/keplertools?branch=main)
[![Documentation Status](https://readthedocs.org/projects/keplertools/badge/?version=latest)](https://keplertools.readthedocs.io/en/latest/?badge=latest)
[![PyPI version](https://badge.fury.io/py/keplertools.svg)](https://badge.fury.io/py/keplertools)
[![Requirements Status](https://requires.io/github/dsavransky/keplertools/requirements.svg?branch=main)](https://requires.io/github/dsavransky/keplertools/requirements/?branch=main)

## Installation

```
pip install SpaceVehiclePropagator
```

To also compile the Cython versions (compiler required, see here for details: https://cython.readthedocs.io/en/latest/src/quickstart/install.html):

```
pip install --no-binary SpaceVehiclePropagator SpaceVehiclePropagator[C]
```

If using a zsh shell (or depending on your specific shell setup), you may need to escape the square brackets (i.e., the last bit of the previous command would be ``SpaceVehiclePropagator\[C\]``.

## Documentation

https://SpaceVehiclePropagator.readthedocs.io/
