# AUBOS

**(*A*xisymmetric *U*nified *B*ackground *O*riented *S*chlieren)**

[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](https://lbesson.mit-license.org/)

The program is designed to analyze background-oriented Schlieren data for axisymmetric objects, primarily using a generalization of the Abel transform. The focus is placed on interpreting data within the unified framework ([Grauer and Steinberg, 2020][GrauerSteinberg20]), thus implementing axisymmetric unified background-oriented Schlieren (AUBOS), and with the use of Bayesian inference and priors. 

### Code contents

The `kernel` package includes functions to generate the typical forward and inverse operators for solving the Abel problem. Inputs vary depending on the operator, with some specifically built for deflectometry measurements, while others apply to the more generic problem and require that the data be transformed prior to use. 

The `tools` package contains miscellaneous functions to aid in analysis. This includes a text-based toolbar function attributed to @sgrauer. 

--------

This code builds on the work by Samuel Grauer (@sgrauer). 

##### References

[Grauer, S. J., & Steinberg, A. M. (2020). Fast and robust volumetric refractive index measurement by unified background-oriented schlieren tomography. Experiments in Fluids, 61(3), 1-17.][GrauerSteinberg20]

[GrauerSteinberg20]: https://link.springer.com/article/10.1007/s00348-020-2912-1
