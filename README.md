---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->



# Bayesian Hierarchical Poisson Models for Multiple Grouped Outcomes with Clustering

<!-- badges: start -->
<!-- badges: end -->

bhpm was developed for the Precision Drug Theraputics: Risk Prediction in Pharmacoepidemiology project as part of a Rutherford Fund Fellowship at Health Data Research (UK), Medical Research Council (UK) award reference MR/S003967/1 (<https://gtr.ukri.org/>).

The goal of bhpm is to investigate associations between multiple outcomes and corresponding patient treatments. bhpm implements Bayesian hierarchical models, which allow a stratification of the population into clusters with similar characteristics,
and which take advantage of known relationships between clinical outcomes, to determine which outcomes are associated with treatments. 

## Installation

You can install the released version of bhpm from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("bhpm")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("rcarragh/bhpm")
```

## Example

This is a basic example which shows how to fit the model:


```r
library(bhpm)
data(demo.cluster.data)
mod.fit <- bhpm.pm(demo.cluster.data, burnin = 100, iter = 200)
#> Memory Model: HIGH
#> MCMC fitting complete.
```
