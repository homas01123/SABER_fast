# S.A.B.E.R (Semi Analytical Bayesian Estimate Retrieval)

A Semi Analytically paramterized aquatic radiative transfer model to primariliy retrieve posterior distributions of Optically Significant Constituents (OSCs), water depth and benthic reflectance from user input of remote sensing reflectance.

Creator and original developer: Soham Mukherjee
Packaging: Raphael Mabit
Maintenance: Soham Mukherjee

## Outline of the mathemtical and physics-based formulation

The user is refered to the publication ``A Semi-Analytical Bayesian Estimate Retrieval (SABER) algorithm for the inversion of Remote Sensing Reflectance in optically deep and shallow waters'' (https://doi.org/10.1002/lom3.70004)

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


## Installation

### Prerequisites

Before installing SABER, ensure you have the following system requirements:

- **R** (>= 3.5.0)
- **pkg-config** 
- **saber-lib** (>= 0.1.0) - C library for fast computations
- **Rtools** (Windows users) or appropriate C++ compiler

### Install from GitHub

Since SABER is currently in development, install directly from GitHub using `devtools`:

```r
# Install devtools if you haven't already
if (!require("devtools")) {
  install.packages("devtools")
}

# Install SABER from GitHub
devtools::install_github("homas01123/SABER_fast")
```

### Install with Dependencies

To install SABER along with suggested packages for full functionality:

```r
# Install with all suggested packages
devtools::install_github("homas01123/SABER_fast", dependencies = TRUE)

# Or install specific suggested packages manually
install.packages(c("knitr", "rmarkdown", "testthat", "ggplot2", "dplyr", "tidyr", "purrr", "tibble", "readr", "magrittr", "BayesianTools"))
```

```

### Verify Installation

Test that SABER is properly installed:

```r
# Load the package
library(SABER)

# Check available forward models
list_forward_models()

# Check available objective functions
list_objective_functions()

# Verify C library integration
build_cache(seq(400, 700, by = 10))
```

### Troubleshooting

**Windows Users:**
- Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) matching your R version
- Ensure Rtools is in your system PATH

**Linux/Mac Users:**
- Install development tools: `sudo apt-get install build-essential` (Ubuntu) or Xcode Command Line Tools (Mac)
- Ensure pkg-config is available: `sudo apt-get install pkg-config`

**Common Issues:**
- If compilation fails, check that all system requirements are met
- For missing saber-lib, refer to the [standalone branch](https://github.com/homas01123/SABER_fast/tree/standalone) for applications

## Outline of the mathemtical and physics-based formulation

The user is refered to the publication ``A Semi-Analytical Bayesian Estimate Retrieval (SABER) algorithm for the inversion of Remote Sensing Reflectance in optically deep and shallow waters'' (htt[...]

## Outline of the code structure

This package follow the recommendations of https://r-pkgs.org/ and the tidyverse style guide https://style.tidyverse.org/.

The code is written with a functional approach. The "business" logic
(the low level functions) are written in file names starting with (`fct_*`),
more generic function are written under (`utils_*`) files. The `fct_*` files usually stores procedures for the low-level computations, that refer to `C` written functions (or its compiled objects)[...]

The central piece of code that stitches together the forward models with the objective function is `objective_factory` high level function.
It allows to easily combine any forward model with any objective functions, even those that a user might add.

For users to add their own forward models and objectives function, we adopted the use of registries.
A registry is an environment in which specific function are stored and can
be retrieved by name in the higher level function arguments.
The three registry currently in use are `.input_preparer_registry`, `.forward_model_registry`, `.objective_function_registry`.

## Running the code

Please refer to the "forward_inverse_basics" vignette for introduction on forward and inverse modeling using SABER.

