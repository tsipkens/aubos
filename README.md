# AUBOS

**(*A*xisymmetric *U*nified *B*ackground-*O*riented *S*chlieren)**

[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](https://lbesson.mit-license.org/)

This program is designed to analyze background-oriented schlieren data for axisymmetric objects, primarily using a generalization of the Abel transform. The focus is placed on interpreting data within the unified framework ([Grauer and Steinberg, 2020][GrauerSteinberg20]), thus implementing axisymmetric unified background-oriented schlieren (AUBOS), and with the use of Bayesian inference and priors. 

## Installation note

This program has a single dependency that are included as submodules: the `cmap` package available at https://github.com/tsipkens/cmap. As a result, this folder will initially be empty. The submodules can be downloaded manually from the above sources and placed in the `cmap/` folder. If cloning using git, clone the repository using 

```shell
git clone git://github.com/tsipkens/aubos --recurse-submodules
```

which will automatically download the submodule. To be used directly, these packages should then be added to the MATLAB path at the beginning of any script using

```Matlab
addpath cmap;
```

Instead of the `cmap` package, one could also replace references in existing scripts to the colormaps that would otherwise be in that package. This would have to include removing any reference to the `cmap_sweep` function (which allows for line plots to sweep through a colormap) that appears in some of the main scripts.   

## Components

This codebase is broken up into a series of packages: 

1. The **transforms** package contain functions explicitly evaluating the Abel and new transform described by Sipkens et al.
2. The **kernel** package includes functions to generate the typical forward and inverse operators for solving the Abel problem. 
3. The **tools** package contains miscellaneous functions to aid in analysis. This includes a text-based toolbar function attributed to @sgrauer. 
4. The **regularization** package contains tools to help during inversion, such as generating prior covariance matrices. 

We refer the reader to individual functions for more information. 

This codebase also contains three classes: 

1. The **Aso** class is used to handle axisymmetric objects that are only defined with respect to radial position (i.e., do not have axial variations). The functions in the +kernel folder are built to evaluate the transforms for these objects. 

2. The **Aso2** class, similarly, is used to handle axisymmetric objects, this time by including radial and axial variations. Solving these problems generally takes much longer than the 1D case considered above. 

3. Finally, the **Camera** class is used to output the ray positions and directions for a pinhole camera, which established a framework by which to expand this representation to include other effects (e.g., lens aberration). 



## Description

The coordinate system used here for the overall axisymmetric schlieren problem is shown below. 

![coord](docs/imgs/01_coordinate.png)

The positive *z*-direction is chosen to proceed forward, away from the camera, and perpindicular to the imaging plane. The origin is placed at the middle of the axisymmetric target object (ASO), such that *z* = 0 represents the distance from the camera lens to the center of the ASO along the imaging axis. 

Projecting axisymmetric objects is typically achieved using the Abel transform, which has a kernel of 

![](https://latex.codecogs.com/svg.latex?\frac{{\partial \delta}}{\partial r} \frac{y_0}{\sqrt{r^2-y_0^2}})

Sipkens et al. (Submitted) derived a new transform, not requiring that the rays passing through the ASO be parallel, which has a kernel of

![](https://latex.codecogs.com/svg.latex?\frac{{\partial \delta}}{\partial r} \frac{y_0}{\sqrt{r^2-y_0^2(1+m_{\text{y}}^2)^{-1}}})

These raw transforms can be evaluated using the functions in the +transforms folder by appending `transform.`  before the function name. For example, the direct, Abel transform can be evaluated using

```Matlab
K = transform.abeld(y0, r_vec);
```

Use of this codebase to evaluate these transforms is demonstrated in the `main_transforms` script. 

### Representing cameras

Imaging inherently requires the use of cameras. Multiple options exist for defining a camera within this program. In any case, one must define the initial trajectories for rays leaving the camera, which will then be used with other components of this codebase. Each ray should be represented by a series of four parameters: 

1. `x0` - The x-position (i.e., the radial position in the imaging plane) at which the ray crosses *z* = 0 (which corresponds to the center of the ASO). 
2. `y0` - Similar to above, this is the y-position (i.e., axial position in the imaging plane) at which the ray crossed *z* = 0.
3. `mx` - The slope of the ray in the *x*-*z* plane. This acts as an input to the transform defined by Sipkens et al. in the associated work. 
4. `my` - The slope of the ray in the *y*-*z* plane. This determines how quickly the ray traverses the axial direction. The `my = 0` case corresponds to rays that do not traverse axially (as would be required for the Abel inversion scenerio).

We provide two examples of how one can define these properties. 

The first involves manually setting the camera properties.  Within the examples provided with this codebase, this is used extensively whenever one wants to focus on the deflection field for only rays in the proximity of the ASO (in other words, ignoring the larger field of view that may be relevant to a real camera). In this case, one can set a camera position and, assuming a pinhole camera, find the trajectory of rays that would original from the pinhole camera and transect the *z* = 0 plane at certain positions. 

While the above treatment is useful within the context of visualizing theoretical deflection fields, as was relevant in generating figures for the associated work by Sipkens et al., more often cameras will be specified instead with a camera origin and focal length. 

--------

### Acknowledgements

This code contains several excerpts from a previous, private codebase by Samuel Grauer (@sgrauer) that were modified for use with this program. 

### References

[Grauer, S. J., & Steinberg, A. M. (2020). Fast and robust volumetric refractive index measurement by unified background-oriented schlieren tomography. Experiments in Fluids, 61(3), 1-17.][GrauerSteinberg20]

[GrauerSteinberg20]: https://link.springer.com/article/10.1007/s00348-020-2912-1
