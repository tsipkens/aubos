# AUBOS

**(*A*xisymmetric *U*nified *B*ackground *O*riented *S*chlieren)**

[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](https://lbesson.mit-license.org/)

This program is designed to analyze background-oriented schlieren data for axisymmetric objects, primarily using a generalization of the Abel transform. The focus is placed on interpreting data within the unified framework ([Grauer and Steinberg, 2020][GrauerSteinberg20]), thus implementing axisymmetric unified background-oriented schlieren (AUBOS), and with the use of Bayesian inference and priors. 

### Installation note

This program has a single dependency that are included as submodules: the `cmap` package available at https://github.com/tsipkens/cmap. As a result, this folder will initially be empty. The submodules can be downloaded manually from the above sources and placed in the `cmap/` folder. If cloning using git, clone the repository using 

```shell
git clone git://github.com/tsipkens/aubos --recurse-submodules
```

which will automatically download the submodule. To be used directly, these packages should then be added to the MATLAB path at the beginning of any script using

```Matlab
addpath('cmap');
```

In the place of the `cmap` package, one could also replace references in existing scripts to the colormaps that would otherwise be in that package.  

# Description
This codebase has a dependency in the form of the `cmap` folder. 

The codebase is broken up into a series of packages: 

1. The `kernel` package includes functions to generate the typical forward and inverse operators for solving the Abel problem. Inputs vary depending on the operator, with some specifically built for deflectometry measurements, while others apply to the more generic problem and require that the data be transformed prior to use. 
2. The `tools` package contains miscellaneous functions to aid in analysis. This includes a text-based toolbar function attributed to @sgrauer. 
3. The `transforms` package contain functions explicitly evaluating the Abel and new transform described by Sipkens et al.
4. The `regularization` package contains tools to help during inversion, such as generating prior covariance matrices. 

We refer the reader to individual functions for more information. 

--------

This code builds on the work by Samuel Grauer (@sgrauer). 

##### References

[Grauer, S. J., & Steinberg, A. M. (2020). Fast and robust volumetric refractive index measurement by unified background-oriented schlieren tomography. Experiments in Fluids, 61(3), 1-17.][GrauerSteinberg20]

[GrauerSteinberg20]: https://link.springer.com/article/10.1007/s00348-020-2912-1
