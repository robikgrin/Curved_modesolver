# Cuved-WG-modesolver
## Overview
Photonic modesolver for curved planar waveguides.

Method based on finite-difference frequency domains (FDFD) method. About FDFD with index averaging you can read in article by [Zhaoming Zhu and Thomas G. Brown](https://www.researchgate.net/publication/24436723_Full-vectorial_finite-difference_analysis_of_microstructure_optical_fibers/citations).

- fully vectorial option with index averaging technique implementation,
- simple structure drawing,
- perfectly-meatched layer **([PML](https://en.wikipedia.org/wiki/Perfectly_matched_layer))** implementation,
- mode profiles visualization (and other plots) via maplotlib.pyplot,
- overlap calculaion

## Examples and features implementations
* [Ex: Mode calculations for curved planar waveguide](#example-mode-calculations-for-curved-planar-waveguide)
* [Feat1: Effective indexes dependency from a width W of a core](#feature-1-effective-indexes-dependency-from-a-width-W-of-a-core)
* [Feat2: Effective indexes dependency from the step of a simulation](#feature-2-effective-indexes-dependency-from-a-step-of-simulation)
* [Feat3: Effective indexes dependency from the distance to the layout edges](#feature-3-effective-indexes-dependency-from-the-distance-to-the-layout-edges)
* [Feat4: Effective indexes dependency from various curvarures of the waveguide](#feature-4-effective-indexes-dependency-from-various-curvatures-of-the-waveguide)


## Introduction
### Main parameters of the planar waveguide and PML
Main parameters of waveguide crossection is on the picture below

<img src="./fiqures/parameters.png " width="600">

Also the picture of parameters of the PML parameters, implemented into our scheme of calculations

<img src="./fiqures/PML_area.png " width="600">

All these parameters ($x_{left}, x_{right}, y_{up}, y_{down}$) needs to be written into into the number of grids, which you need to pave, f.e. $x_{left} = 100$ means, that the width of the PML at the left will be $100 \times d \xi$ (read about $d \xi$ below).

### Other parameters of simulation
- wavelength $\lambda$ (``lambda``)
- steps for finite-differences in each direction: $d \xi$ (``d_xi``) in horizontal direction, $d \eta$ (``d_eta``) - in vertical direction
- curvature of waveguide $\kappa$ (``kappa``)

## Example: Mode calculations for straight/curved planar waveguide
Let's see how programm will work with default settings. The default setting is:

* wavelength: 1.55 $\mu m$
* clad refractive index: 1.4444 ($SiO_2$ cladding)
* core refractive index: 3.4755 ($Si$ core)
* grid step size in $\xi$ direction: 0.02 $\mu m$
* grid step size in $\eta$ direction: 0.02 $\mu m$
* width of $Si$ core: 2 $\mu m$
* height of $Si$ core: 0.22 $\mu m$
* distance to the left border of simulation: 2 $\mu m$
* distance to the right border of simulation: 2 $\mu m$
* distance to the upper border of simulation: 2 $\mu m$
* distance to the down border of simulation: 2 $\mu m$
* curvature value: 0 $\mu m^{-1}$

In our script we calculating first 2 modes of waveguide.

### Python script

```python
import curved_modesolver.py as cms
import numpy as np

keke = rect_WG() #constructing the object with default parameters


NPML_l = [100, 100, 100, 100]

keke.FDE(2, NPML_l)

#Visualization
keke.draw_structure()
keke.draw_field(1, 'norm')
keke.draw_field(2, 'norm')
```

#### Structure

<img src="./fiqures/ex_def_structure.png " width="600">

#### Modes

<img src="./fiqures/ex_def_1_mode.png " width="600">

<img src="./fiqures/ex_def_2_mode.png " width="600">

## Feature 1: Effective indexes dependency from a width W of a core

bla-bla-bla

## Feature 2: Effective indexes dependency from the step of a simulation

bla-bla-bla

## Feature 3: Effective indexes dependency from the distance to the layout edges

bla-bla-bla

## Feature 4: Effective indexes dependency from various curvarures of the waveguide

bla-bla-bla


## Future features
- :star:Any possible crossection profile calculations

And maybe something more, we will see :smirk: 