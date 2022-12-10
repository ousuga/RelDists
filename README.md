
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![R-CMD-check](https://github.com/ousuga/RelDists/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ousuga/RelDists/actions/workflows/R-CMD-check.yaml)
[![Travis build
status](https://travis-ci.org/ousuga/RelDists.svg?branch=master)](https://travis-ci.org/ousuga/RelDists)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/ousuga/RelDists?branch=master&svg=true)](https://ci.appveyor.com/project/ousuga/RelDists)
[![Lifecycle:
maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html#maturing)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version-ago/RelDists)](https://cran.r-project.org/package=RelDists)
[![CRAN_Download_Badge](http://cranlogs.r-pkg.org/badges/RelDists)](https://cran.r-project.org/package=RelDists)
<!-- badges: end -->

# RelDists <img src="man/figures/RelDists4.3_gris.png" align="right" height="200" align="right"/>

In this package are available multiple useful distributions for
reliability analysis. With this package it is possible to estimate
parameters and fit regression models within GAMLSS framework.

## Installation

To install the `RelDists` package you need to install `devtools` package
before, follow the instructions below:

``` r
if (!require('devtools')) install.packages('devtools')
devtools::install_github('ousuga/RelDists', force=TRUE)
require(RelDists)
```

You can visit the [package website](https://ousuga.github.io/RelDists/)
to explore the vignettes (articles) and function reference.
