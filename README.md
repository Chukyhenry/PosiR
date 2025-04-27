
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PosiR <img src="man/figures/logo.png" align="right" height="140"/>

[![R-CMD-check](https://github.com/Chukyhenry/PosiR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Chukyhenry/PosiR/actions)

**PosiR** provides valid post-selection inference (PoSI) for linear
regression models, based on simultaneous confidence intervals using
bootstrap-based max-t statistics.

It implements Algorithm 1 from [Kuchibhotla et
al. (2022)](http://dx.doi.org/10.1146/annurev-statistics-100421-044639),
designed to give valid coverage even after variable selection.

------------------------------------------------------------------------

## Installation

``` r
# Install development version from GitHub
devtools::install_github("Chukyhenry/PosiR")
#> Using GitHub PAT from the git credential store.
#> Downloading GitHub repo Chukyhenry/PosiR@HEAD
#> ── R CMD build ──────────────────────────────────────────
#>      checking for file ‘/private/var/folders/4f/jz1fz83n0tg9csh3mc718cs80000gn/T/RtmpnZ82nH/remotes1362b55cc26a4/Chukyhenry-PosiR-ff3be1e/DESCRIPTION’ ...  ✔  checking for file ‘/private/var/folders/4f/jz1fz83n0tg9csh3mc718cs80000gn/T/RtmpnZ82nH/remotes1362b55cc26a4/Chukyhenry-PosiR-ff3be1e/DESCRIPTION’
#>   ─  preparing ‘PosiR’:
#>    checking DESCRIPTION meta-information ...  ✔  checking DESCRIPTION meta-information
#>   ─  checking for LF line-endings in source and make files and shell scripts
#>   ─  checking for empty or unneeded directories
#>    Omitted ‘LazyData’ from DESCRIPTION
#>   ─  building ‘PosiR_0.1.0.tar.gz’
#>      
#> 
#> Installing package into '/Users/henrychukwuma/Library/Caches/org.R-project.R/R/renv/library/PosiR-741335c4/macos/R-4.4/aarch64-apple-darwin20'
#> (as 'lib' is unspecified)
#> Warning in i.p(...): installation of package
#> '/var/folders/4f/jz1fz83n0tg9csh3mc718cs80000gn/T//RtmpnZ82nH/file1362b4dfd04fd/PosiR_0.1.0.tar.gz'
#> had non-zero exit status
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
devtools::load_all()
#> ℹ Loading PosiR
# Simulate some data
set.seed(123)
X <- matrix(rnorm(100 * 3), 100, 3)
colnames(X) <- c("X1", "X2", "X3")
y <- 1 + X[,1] * 0.5 + rnorm(100)

# Define models (column indices)
Q <- list(model1 = 1:2, model2 = 1:4)

# Run PoSI
result <- simultaneous_ci(X, y, Q, B = 500, verbose = FALSE)

# View results
print(result$intervals)
#>   model_id coefficient_name    estimate      lower
#> 1   model1      (Intercept)  0.96831201  0.7019780
#> 2   model1               X1  0.44983825  0.1861455
#> 3   model2      (Intercept)  0.98067390  0.7074605
#> 4   model2               X1  0.44454938  0.1809392
#> 5   model2               X2  0.04621817 -0.2186037
#> 6   model2               X3 -0.05738700 -0.3736027
#>       upper psi_hat_nqj    se_nqj
#> 1 1.2346460    1.084196 0.1041248
#> 2 0.7135310    1.062799 0.1030922
#> 3 1.2538873    1.140929 0.1068143
#> 4 0.7081596    1.062134 0.1030599
#> 5 0.3110400    1.071920 0.1035336
#> 6 0.2588287    1.528346 0.1236263
plot(result)
```

<img src="man/figures/README-example-1.png" width="100%" />

## Learn More

Vignette: **vignette(“introduction_to_posir”)** [Online
Documentation](https://chukyhenry.github.io/PosiR/) Source Paper:
**Kuchibhotla et al. (2022)**
