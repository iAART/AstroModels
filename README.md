[![2211.07469](https://img.shields.io/badge/arXiv-2411.17884-b31b1b.svg)](https://arxiv.org/abs/2411.17884) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/iAART/aart/License.txt) [![GitHub repo stars](https://img.shields.io/github/stars/iAART/AstroModels?style=social)](https://github.com/iAART/AstroModels)
[![DOI](https://zenodo.org/badge/873706198.svg)](https://doi.org/10.5281/zenodo.14629535)


# Multi-frequency Models for AART
_______
Desire, T., Cardenas-Avendano, A., & Chael, A. 2024, arXiv e-prints, [arXiv:2411.17884](http://arxiv.org/abs/2411.17884)
_______
This repository is an extended version of AART, a full detailed description of which is at https://github.com/iAART/aart.

These modifications focus on various physical parameters affecting the AART emission profile behavior. The simulation includes calculations of temperature, density, magnetic fields, and specific intensity, tailored for astrophysical contexts.

_______

We also request that modifications or extensions leading to a scientific publication be made public as free software. 

<center> <em>Feel free to use images and movies produced with this code (with attribution) for your next presentation! </em> </center>

_______
![GitHub last commit](https://img.shields.io/github/last-commit/iAART/AstroModels)
_______


## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Key Features](#key-features)
- [Parameters](#parameters)
- [Example](#examples)

_______
### Installation

To use this simulation, clone the repository and install the required packages:

```bash
git clone https://github.com/iAART/AstroModels.git
cd AstroModels
pip install -r requirements.txt
```

Ensure you have Python 3.x and the following libraries installed:

    NumPy
    SciPy
    Astropy
    lmfit
_______
### Key Features

   * **Radial Profiles:**
   All major additions are contained in AstroModels/aart_func/rprofs_f.py
    
   * Temperature, Density Profiles, and Magnetic Field Strength: Calculate the electron temperature and density as functions of distance from the black hole, as well as the Magnetic field strength.
   * Synchrotron Emission: Models synchrotron radiation for different configurations and astrophysical parameters.
   * Noisy Density and Temperature Profiles: Incorporates noise in density and temperature for realistic simulations.
_______
### Parameters

The following parameters can be adjusted to modify the simulation:

|Parameter| Description                                                       |    Base Value |
| ------- | ----------------------------------------------------------------- | ------------- |
| nu      | Observation frequency                                             |          230e9|
| mass    | blackhole mass                                                    |     6.5e9 Msun| 
| scaleh  | plasma disk scale height                                          |            0.5|
| theta_b | angle between magnetic field and wave vector, ignored if not fixed|     60 degrees|
| rb0     | distance at which power laws take base value                      |           5r_g|
| nth0    | Base value for the density power law                              | 3.99e5 1/cm^-3|
| te0     | Base value for temperature power law                              |          3e10K|
| b0      | Base value for the magnetic field power law                       |        8 Gauss|
| pdens   | Sets decay speed for the density power law                        |           -0.7|
| ptemp   | Sets decay speed for the temperature power law                    |             -1|
| pmag    | Sets decay speed for the magnetic field power law                 |           -1.5|
| nscale  | inosiy perturbation scale                                         |             .3|
   
* Function keys: The choice to add Inoisy purturbations to the density, temperature, and/or magnetic field strength power laws.

Refer to the code comments for detailed descriptions of each parameter and its units.

_______
### Examples

And example run of the code can be found as "ExampleMultiFrequency.ipynb"
