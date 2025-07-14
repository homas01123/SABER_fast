# S.A.B.E.R (Semi Analytical Bayesian Estimate Retrieval)

A Semi Analytically paramterized aquatic radiative transfer model to primariliy retrieve posterior distributions of Optically Significant Constituents (OSCs), water depth and benthic reflectance from user input of remote sensing reflectance.

Creator and original developer: Soham Mukherjee
Packaging: Raphael Mabit
Maintenance: Soham Mukherjee

## Outline of the mathemtical and physics-based formulation

The user is refered to the publication ``A Semi-Analytical Bayesian Estimate Retrieval (SABER) algorithm for the inversion of Remote Sensing Reflectance in optically deep and shallow waters'' (https:/)

## Outline of the code structure

This package follow the recommendations of https://r-pkgs.org/ and the tidyverse style guide https://style.tidyverse.org/.

The code is written with a functional approach. The "business" logic
(the low level functions) are written in file names starting with (`fct_*`),
more generic function are written under (`utils_*`) files. The `fct_*` files usually stores procedures for the low-level computations, that refer to `C` written functions (or its compiled objects) inside the `src`.

The central piece of code that stitches together the forward models with the objective function is `objective_factory` high level function.
It allows to easily combine any forward model with any objective functions, even those that a user might add.

For users to add their own forward models and objectives function, we adopted the use of registries.
A registry is an environment in which specific function are stored and can
be retrieved by name in the higher level function arguments.
The three registry currently in use are `.input_preparer_registry`, `.forward_model_registry`, `.objective_function_registry`.

## Running the code
