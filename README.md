# Multi-frequency Models for AART
This repository is an extended version of AART, a full detailed description of which is at https://github.com/iAART/aart.

These modifications focus on various physical parameters affecting the AART emission profile behavior. The simulation includes calculations of temperature, density, magnetic fields, and specific intensity, tailored for astrophysical contexts.

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Key Features](#key-features)
- [Parameters](#parameters)
- [Example](#examples)

## Installation

To use this simulation, clone the repository and install the required packages:

```bash
git clone https://github.com/yourusername/black-hole-simulation.git
cd black-hole-simulation
pip install -r requirements.txt
```

Ensure you have Python 3.x and the following libraries installed:

    NumPy
    SciPy
    Astropy
    lmfit

Key Features

    Temperature, Density Profiles, and Magnetic Field Strength: Calculate the electron temperature and density as functions of distance from the black hole, as well as the Magnetic field strength.
    Synchrotron Emission: Models synchrotron radiation for different configurations and astrophysical parameters.
    Noisy Density and Temperature Profiles: Incorporates noise in density and temperature for realistic simulations.

Parameters

The following parameters can be adjusted to modify the simulation:

    Mass: Mass of the black hole.
    Observation Frequency: Frequency of observation in Hz.
    Scale Height: The slope of the accretion disk height vs radius.
    Initial Conditions: Various coefficients and exponents related to density, temperature, and magnetic field strength power laws.
    Function keys: The choice to add Inoisy purturbations to the density, temperature, and/or magnetic field strength power laws.

Refer to the code comments for detailed descriptions of each parameter and its units.
Examples

And example run of the code can be found as "ExampleMultiFrequency.ipynb"
