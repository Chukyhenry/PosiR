
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PosiR <img src="man/figures/logo.png" align="right" height="140"/>

[![R-CMD-check](https://github.com/Chukyhenry/PosiR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Chukyhenry/PosiR/actions)

**PosiR** provides tools for post-selection inference (PoSI) in linear
regression models. PoSI addresses the challenge of performing valid
statistical inference after model selection, ensuring that confidence
intervals maintain their nominal coverage probability (e.g., 95%) even
when the model is chosen based on the data. The package implements
simultaneous confidence intervals using bootstrap-based max-t
statistics, following Algorithm 1 from Kuchibhotla et al. (2022).

## Installation

You can install the development version of `PosiR` from GitHub:

``` r
# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
# Install PosiR
devtools::install()
#> ── R CMD build ───────────────────────────────────────────────────
#>      checking for file ‘/Users/henrychukwuma/Documents/Spring 2025/STAT 8054/Project/PosiR/DESCRIPTION’ ...  ✔  checking for file ‘/Users/henrychukwuma/Documents/Spring 2025/STAT 8054/Project/PosiR/DESCRIPTION’
#>   ─  preparing ‘PosiR’:
#>    checking DESCRIPTION meta-information ...  ✔  checking DESCRIPTION meta-information
#>   ─  checking for LF line-endings in source and make files and shell scripts
#>   ─  checking for empty or unneeded directories
#>    Removed empty directory ‘PosiR/man/figures’
#> ─  building ‘PosiR_0.1.0.tar.gz’
#>      
#> Running /Library/Frameworks/R.framework/Resources/bin/R CMD \
#>   INSTALL \
#>   /var/folders/4f/jz1fz83n0tg9csh3mc718cs80000gn/T//Rtmp2lOOaq/PosiR_0.1.0.tar.gz \
#>   --install-tests 
#> * installing to library ‘/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library’
#> * installing *source* package ‘PosiR’ ...
#> ** using staged installation
#> ** R
#> ** tests
#> ** byte-compile and prepare package for lazy loading
#> ** help
#> *** installing help indices
#> ** building package indices
#> ** installing vignettes
#> ** testing if installed package can be loaded from temporary location
#> ** testing if installed package can be loaded from final location
#> ** testing if installed package keeps a record of temporary installation path
#> * DONE (PosiR)

# Optional dependencies for vignette and examples
install.packages(c("dplyr", "pbapply"))
#> 
#> The downloaded binary packages are in
#>  /var/folders/4f/jz1fz83n0tg9csh3mc718cs80000gn/T//Rtmp2lOOaq/downloaded_packages
```

## Example: Simultaneous Confidence Intervals

This example demonstrates how to use `simultaneous_ci()` to compute
simultaneous confidence intervals for regression coefficients across a
set of models:

``` r
library(PosiR)

# Simulate data
set.seed(123)
X <- matrix(rnorm(100 * 3), 100, 3)
colnames(X) <- c("X1", "X2", "X3")
y <- 1 + X[, "X1"] * 0.5 + rnorm(100)  # True intercept = 1, X1 coefficient = 0.5

# Define model universe (column indices of X)
Q <- list(
  model1 = 1:2,  # Model with X1, X2
  model2 = 1:3   # Model with X1, X2, X3
)

# Compute simultaneous confidence intervals
result <- simultaneous_ci(X, y, Q, B = 500, verbose = FALSE)

# View results
print(result$intervals)
#>   model_id coefficient_name   estimate      lower     upper
#> 1   model1      (Intercept) 0.96831201  0.7198033 1.2168207
#> 2   model1               X1 0.44983825  0.2037940 0.6958825
#> 3   model2      (Intercept) 0.97292290  0.7230406 1.2228052
#> 4   model2               X1 0.45219170  0.2012421 0.7031413
#> 5   model2               X2 0.04485171 -0.1971332 0.2868366
#>   psi_hat_nqj    se_nqj
#> 1    1.084196 0.1041248
#> 2    1.062799 0.1030922
#> 3    1.096215 0.1047003
#> 4    1.105600 0.1051475
#> 5    1.028019 0.1013913

# Plot the intervals
plot(result, main = "Simultaneous Confidence Intervals", las.labels = 1)
```

<img src="man/figures/README-example-1.png" width="100%" /> \##
Interpretation

The output `result$intervals` provides the coefficient estimates and
simultaneous 95% confidence intervals for each model in `Q`. For
example:

The `(Intercept)` and `X1` intervals in `model1` should contain their
true values (1 and 0.5, respectively).

The intervals are wider than naive intervals to account for model
selection uncertainty, ensuring valid coverage across all models in `Q`.

## Learn More

Vignette: **vignette(“Vignette”)** Source Paper: Kuchibhotla, A.,
Kolassa, J., & Kuffner, T. (2022). Post-selection inference. Annual
Review of Statistics and Its Application, 9(1), 505–527. DOI:
10.1146/annurev-statistics-100421-044639.
