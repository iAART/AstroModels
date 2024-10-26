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

    |Parameter| Description                                                       |    Base Value |
    |---------|-------------------------------------------------------------------|---------------|
    | nu      | Observation frequency                                             |          230e9|
    | mass    | blackhole mass                                                    |     6.5e9 Msun| 
    | scaleh  | plasma disk scale height                                          |            0.5|
    | theta_b | angle between magnetic field and wave vector, ignored if not fixed|     60 degrees|
    | rb0     | distance at which power laws take base value                      |            5rg|
    | nth0    | Base value for the density power law                              | 3.99e5 1/cm^-3|
    | te0     | Base value for temperature power law                              |          3e10K|
    | b0      | Base value for the magnetic field power law                       |        8 Gauss|
    | pdens   | Sets decay speed for the density power law                        |              1|
    | ptemp   | Sets decay speed for the temperature power law                    |              1|
    | pmag    | Sets decay speed for the magnetic field power law                 |            1.5|
    | nscale  | inosiy perturbation scale                                         |             .3|
    
    Function keys: The choice to add Inoisy purturbations to the density, temperature, and/or magnetic field strength power laws.

Refer to the code comments for detailed descriptions of each parameter and its units.
Examples

And example run of the code can be found as "ExampleMultiFrequency.ipynb"
