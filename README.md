# Cuved-WG-modesolver
## Overview
Photonic modesolver for curved planar waveguides.

Method based on finite-difference frequency domains (FDFD) method. About FDFD you can read in article by [Zhaoming Zhu and Thomas G. Brown](https://www.researchgate.net/publication/24436723_Full-vectorial_finite-difference_analysis_of_microstructure_optical_fibers/citations). 

- fully vectorial option with index averaging technique implementation,
- simple structure drawing,
- perfectly-meatched layer **(PML)** implementation,
- mode profiles visualization (and other plots) via maplotlib.pyplot,
- overlap calculaion
## Feature features
- [ ] Any possible crossection profile calculations

And maybe something more, we will see. 

## Introduction
### Main parameters of the planar waveguide and PML
Main parameters of waveguide crossection is on the picture below
<img src="./fiqures/parameters.png " width="600">

Also the picture of parameters of the PML parameters, implemented into our scheme of calculations
<img src="./fiqures/PML_area.png " width="600">

All these parameters needs to be written into into the number of grids, which you need to pave, f.e. $x_{left} = 100$ means, that the width of the PML at the left will be $100 \times d \xi$ (read about $d \xi$ below).

### Other parameters of simulation
- wavelength ``lambda``
- steps for finite-differences in each direction: ``d_xi`` in horizontal direction, ``d_eta`` - in vertical direction
- curvature of waveguide ``kappa``

## Examples and features implementations
* [Ex: Mode calculations for curved planar waveguide](#example-mode-calculations-for-curved-planar-waveguide)
* [Feat1: Effective indexes dependency from a width W of a core](#feature-1-effective-indexes-dependency-from-a-width-W-of-a-core)
* [Feat2: Effective indexes dependency from the step of a simulation](#feature-2-effective-indexes-dependency-from-a-step-of-simulation)
* [Feat3: Effective indexes dependency from the distance to the layout edges](#feature-3-effective-indexes-dependency-from-the-distance-to-the-layout-edges)
* [Feat4: Effective indexes dependency from various curvarures of the waveguide](#feature-4-effective-indexes-dependency-from-various-curvatures-of-the-waveguide)

## Example: Mode calculations for curved planar waveguide

bla-bla-bla

## Feature 1: Effective indexes dependency from a width W of a core

bla-bla-bla

## Feature 2: Effective indexes dependency from the step of a simulation

bla-bla-bla

## Feature 3: Effective indexes dependency from the distance to the layout edges

bla-bla-bla

## Feature 4: Effective indexes dependency from various curvarures of the waveguide

bla-bla-bla
