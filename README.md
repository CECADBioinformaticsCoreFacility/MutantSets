
# MutantSets

<!-- badges: start -->
[![Codecov test coverage](https://codecov.io/gh/RichardJActon/MutantSets/branch/master/graph/badge.svg)](https://codecov.io/gh/RichardJActon/MutantSets?branch=master)
[![Travis build status](https://travis-ci.com/RichardJActon/MutantSets.svg?branch=master)](https://travis-ci.com/RichardJActon/MutantSets)
<!-- badges: end -->

The goal of MutantSets is to permit the exploration of the results of whole genome sequencing and mutation calling in *C. elegans*, with the goal of identify candidate mutations responsible for phenotypes in genetic screens though mapping by sequencing.

## Installation

You can install the released version of MutantSets from [github](https://github.com/RichardJActon/MutantSets) with:

``` r
remotes::install_github("RichardJActon/MutantSets")
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
