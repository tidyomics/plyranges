
<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Build Status](https://travis-ci.org/sa-lee/plyranges.svg?branch=master)](https://travis-ci.org/sa-lee/plyranges) [![Coverage Status](https://img.shields.io/codecov/c/github/sa-lee/plyranges/master.svg)](https://codecov.io/github/sa-lee/plyranges?branch=master)

plyranges
=========

Introduction
------------

The `plyranges` package provides a [fluent](https://en.wikipedia.org/wiki/Fluent_interface) interface for performing common genomic data analysis tasks. One goal of this package is to decrease the learning curve for Bioconductor and S4 classes (especially for new users) by providing a consistent interface to common classes, *IRanges*, *GRanges* and *SummarizedExperiment* (not yet implemented).

To achieve this we have constructed a grammar of genomic data manipulation based on `dplyr` and the Bioconductor packages that define the aformentioned classes.

All of the methods defined in `plyranges` are type-consistent and chainable. As a consequence, users familiar with `dplyr` should be able to read and understand `plyranges` code without difficulty.

*Note this package is under heavy development*

Installation
============

Currently only the development version is available:

``` r
source("https://bioconductor.org/biocLite.R")
biocLite("sa-lee/plyranges")

# alternatively with the devtools package
library(devtools)
install_github("sa-lee/plyranges")
```

Getting started
===============

See the [introduction vignette](http://stuartlee.org/plyranges/articles/01-introduction-plyranges.html) for an overview of `plyranges` features. For mapping between `plyranges`, `bedtools` and `GenomicRanges` see the [data wrangling vignette](http://stuartlee.org/plyranges/articles/02-bedtools-examples.html).
