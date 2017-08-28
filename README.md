<!-- README.md is generated from README.Rmd. Please edit that file -->
plyranges
=========

The `plyranges` [1] package provides a (fluent)\[<https://en.wikipedia.org/wiki/Fluent_interface>\] interface for performing common genomic data analysis tasks. One goal of this package is to decrease the learning curve for Bioconductor and S4 classes (especially for new users) by providing a consistent interface to common classes, *IRanges*, *GRanges* and *SummarizedExperiment*.

To achieve this we have constructed a grammar of genomic data manipulation based on `dplyr` and the Bioconductor packages that define the aformentioned classes. All the methods defined in `plyranges` are type-consistent and chainable. As a consequence, users familiar with `dplyr` should be able to read and understand `plyranges` code without difficulty.

*Note this package is under heavy development*

[1] tentative title
