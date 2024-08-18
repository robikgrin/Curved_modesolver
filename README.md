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
* [Ex1: Mode calculations for straight planar waveguide with default parameters](#example-1-mode-calculations-for-straight-planar-waveguide-with-default-parameters)
* [Ex2: Mode calculations for curved planar waveguide](#example-2-mode-calculations-for-curved-planar-waveguide)
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

## Example 1: Mode calculations for straight planar waveguide with default parameters
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

keke = rect_WG() #constructing the object with default parameters


NPML_l = [100, 100, 100, 100]

keke.FDE(2, NPML_l)

#Visualization
keke.draw_structure()
keke.draw_field(1, 'norm')
keke.draw_field(2, 'norm')
```

#### Structure

<img src="./fiqures/ex1/ex1_def_structure.png " width="600">

#### Modes

<img src="./fiqures/ex1/ex1_def_1_mode.png " width="800">

<img src="./fiqures/ex1/ex1_def_2_mode.png " width="800">

## Example 2: Mode calculations for curved planar waveguide
In this script we also will calculate 2 modes. The settings of simulation:

* wavelength: 1.55 $\mu m$
* clad refractive index: 1.4444 ($SiO_2$ cladding)
* core refractive index: 3.4755 ($Si$ core)
* grid step size in $\xi$ direction: 0.02 $\mu m$
* grid step size in $\eta$ direction: 0.02 $\mu m$
* width of $Si$ core: 2 $\mu m$
* height of $Si$ core: 0.22 $\mu m$
* distance to the left border of simulation: 5 $\mu m$
* distance to the right border of simulation: 3 $\mu m$
* distance to the upper border of simulation: 3 $\mu m$
* distance to the down border of simulation: 3 $\mu m$
* curvature value:  $\mu m^{-1}$

```python
import curved_modesolver.py as cms

#Parameters initialization
wavelength = 1.55E-6
n_clad = 1.4444
n_core = 3.4755
d_xi = 2E-8
d_eta = 2E-8
W = 2E-6
H = 2.2E-7
delta_l = 5E-6
delta_r = 5E-6
delta_u = 3E-6
delta_d = 3E-6
kappa = 0.1E6

obj = rect_WG(wavelength, n_clad, n_core, d_xi, d_eta, W, H, delta_l, delta_r, delta_u, delta_d, kappa)

NPML_l = [100, 100, 100, 100]
obj.FDE(2, NPML_l)

#Visualization
obj.draw_structure()
obj.draw_field(1, 'norm')
obj.draw_field(1, 'log')
obj.draw_field(2, 'norm')
obj.draw_field(2, 'log')
```

#### Structure

<img src="./fiqures/ex2/ex2_def_structure.png " width="600">

#### Modes
Norm scale
<img src="./fiqures/ex2/ex2_def_1_mode_norm.png " width="800">

Log scale
<img src="./fiqures/ex2/ex2_def_1_mode_log.png " width="800">

Norm scale
<img src="./fiqures/ex2/ex2_def_2_mode_norm.png " width="800">

Log scale
<img src="./fiqures/ex2/ex2_def_2_mode_log.png " width="800">

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