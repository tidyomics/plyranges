---
name: Bug report or feature request
about: Describe a bug you've seen or make a case for a new feature
title: "[BUG] Your bug or feature request"
labels: ''
assignees: ''
---

Please briefly describe your problem and what output you expect. If you have a question, please don't use this form. Instead, ask on <https://support.bioconductor.org/> using the appropriate tag(s) including one for this package.

## Context

Provide some context for your bug report or feature request. This could be the:

* link to raw code, example: https://github.com/lcolladotor/osca_LIIGH_UNAM_2020/blob/master/00-template.Rmd#L24-L28
* link to a commit, example: https://github.com/lcolladotor/osca_LIIGH_UNAM_2020/commit/6aa30b22eda614d932c12997ba611ba582c435d7
* link to a line of code inside a commit, example: https://github.com/lcolladotor/osca_LIIGH_UNAM_2020/commit/6aa30b22eda614d932c12997ba611ba582c435d7#diff-e265269fe4f17929940e81341b92b116R17
* link to code from an R package, example: https://github.com/LieberInstitute/spatialLIBD/blob/master/R/run_app.R#L51-L55

## Code

Include the code you ran and comments

```R
## prompt an error
stop('hola')

## check the error trace
traceback()
```

## Small reproducible example

If you copy the lines of code that lead to your error, you can then run [`reprex::reprex()`](https://reprex.tidyverse.org/reference/reprex.html) which will create a small website with code you can then easily copy-paste here in a way that will be easy to work with later on.

```R
## prompt an error
stop('hola')
#> Error in eval(expr, envir, enclos): hola

## check the error trace
traceback()
#> No traceback available
```


## R session information

Remember to include your full R session information.

```R
options(width = 120)
sessioninfo::session_info()
```

The output of `sessioninfo::session_info()` includes relevant GitHub installation information and other details that are missed by `sessionInfo()`.
