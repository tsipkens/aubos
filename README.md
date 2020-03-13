# Axis-symmetric, unified background-oriented Schlieren (AUBOS) 

The program is designed to analyze background-oriented Schlieren data for axis-symmetrical objects using the Abel transform and related operations. The focus is placed on interpreting data within the unified framework (Grauer and Steinberg, 2020][GrauerSteinberg20]), thus implementing axis-symmetric unified background-oriented Schlieren (AUBOS). 

### Code contents

The `abel` package includes functions to generate the typical forward and inverse operators for solving the Abel problem. Inputs vary depending on the operator, with some specifically built for deflectometry measurements, while others apply to the more generic problem and require that the data be transformed prior to use. 

The `tools` package contains miscellaneous functions to aid in analysis. This includes a text-based toolbar function attributed to @sgrauer. 

--------

This code builds on the work by Samueal Grauer (@sgrauer). 

##### References

[Grauer, S. J., & Steinberg, A. M. (2020). Fast and robust volumetric refractive index measurement by unified background-oriented schlieren tomography. Experiments in Fluids, 61(3), 1-17.][GrauerSteinberg20]

[GrauerSteinberg20]: https://link.springer.com/article/10.1007/s00348-020-2912-1
