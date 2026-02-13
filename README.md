# Gamma Generalized Inference Project

This repository contains implementation and examples for generalized
p-values and confidence intervals for Gamma and log-Gamma quantiles
using Monte Carlo methods based on:

Weerahandi, S. & Gamage, J. (2016).\
*A general method of inference for two-parameter continuous
distributions*.\
Communications in Statistics---Theory and Methods, 45(9), 2612--2625.\
DOI: 10.1080/03610926.2014.887109.

------------------------------------------------------------------------

## Files in this Repository

-   **gamma_generalized_pvalues.Rmd**\
    R Markdown source file containing:
    -   Method description
    -   Implementation of Monte Carlo T-ratio CDF
    -   Generalized p-value for log-Gamma quantiles
    -   Confidence interval construction
    -   Empirical coverage simulation example
-   **gamma_generalized_pvalues.html**\
    Knitted HTML output generated from the R Markdown file.

------------------------------------------------------------------------

## Main Function

The primary function of interest is:

``` r
CILogGamQuan(data, p = 0.9, cover = 0.95, Nmc = 500, Ncdf = 500)
```

It computes a generalized confidence interval for a log-Gamma quantile.

------------------------------------------------------------------------

## Reproducing Results

1.  Open `gamma_generalized_pvalues.Rmd` in RStudio.
2.  Click **Knit** to generate the HTML file.
3.  Adjust Monte Carlo parameters (`Nmc`, `Ncdf`, and simulation
    repetitions `R`) to trade off between speed and accuracy.

------------------------------------------------------------------------

## Notes

-   Computation can be intensive due to repeated `uniroot()` calls.
-   Increase Monte Carlo sizes for more stable empirical coverage
    results.
-   Smaller values are recommended when testing or knitting the
    document.

## Usage 

set.seed(123)
x <- rgamma(20, shape = 2, scale = 5)
CILogGamQuan(x, p = 0.9, cover = 0.95, Nmc = 500, Ncdf = 500)

