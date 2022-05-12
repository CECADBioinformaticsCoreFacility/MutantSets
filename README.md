
# MutantSets

<!-- badges: start -->
[![Codecov test coverage](https://codecov.io/gh/CECADBioinformaticsCoreFacility/MutantSets/branch/master/graph/badge.svg)](https://codecov.io/gh/CECADBioinformaticsCoreFacility/MutantSets?branch=master)
[![Travis build status](https://travis-ci.com/CECADBioinformaticsCoreFacility/MutantSets.svg?branch=master)](https://travis-ci.com/CECADBioinformaticsCoreFacility/MutantSets)
<!-- badges: end -->

The goal of MutantSets is to permit the exploration of the results of whole genome sequencing and mutation calling in *C. elegans*, with the goal of identify candidate mutations responsible for phenotypes in genetic screens though mapping by sequencing.

## Installation

You can install the released version of MutantSets from [github](https://github.com/CECADBioinformaticsCoreFacility/MutantSets) with:

``` r
remotes::install_github("CECADBioinformaticsCoreFacility/MutantSets")
```

NB on more recent R versions there may be issues install the `vcfR` package if this is an issue when you attempt to install this package please try installing `vcfR` from source with:

```r
remotes::install_github("knausb/vcfR")
```

## Example

To start the app run:

``` r
MutantSets::launchApp()
```

For instructions on how to prepare data for use in this app see the vignette:

```r
vignette("Generating-input", package = "MutantSets")
```

## Run on Renku

https://renkulab.io/projects/racton/mutantsets

You can Run MutantSets without installing it on your local machine by going to it's page on Renkulab.io and hitting start, it will take a minute to get spun up but should work just fine after that. Click open to view in full screen.
