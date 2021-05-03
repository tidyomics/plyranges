# this enables the use of `::` so we can dispatch on the
# right generics

#' @importFrom methods getGeneric getGenerics
bioc_generics <- function(pkgs = c("BiocGenerics", "IRanges", "S4Vectors")) {
  pkgs <- lapply(pkgs, asNamespace)
  generics <-  getGenerics(pkgs)
  fn <- mapply(getGeneric,
               f = generics@.Data,
               package = generics@package,
               SIMPLIFY = FALSE)
  Filter(function(x) {
    fun = try(x@default, silent = TRUE)
    if (is(fun, "try-error")) FALSE
    !is.primitive(fun)
  },
  fn)
}

# bind the generics to a child env, so a quosure has access to them

scope_plyranges <- function(env, generics) {
  tail <- env
  nms <- character(0)
  # recurse through parent environments until we get to the empty env
  while (!identical(tail, rlang::empty_env())) {
    env_nms <- rlang::env_names(tail)
    nms <- unique(c(nms, intersect(names(generics), env_nms)))
    tail <- rlang::env_parent(tail)
  }
  nms <- setdiff(nms, rlang::env_names(rlang::global_env()))
  generics <- generics[nms]

  child <- rlang::child_env(env, !!!generics)
  child
}

# Given a set of quosures bind the generics as a child environment
quos_set_env <- function(quos, generics = bioc_generics()) {
  for (i in seq_along(quos)) {
    quos[[i]] <- rlang::quo_set_env(
      quos[[i]],
      env = scope_plyranges(rlang::quo_get_env(quos[[i]]), generics)
    )
  }
  quos
}


set_dots_unnamed <- function(...) {
  dots <- rlang::enquos(...)
  quos_set_env(dots)
}

set_dots_named <- function(...) {
  dots <- rlang::enquos(..., .named = TRUE)
  quos_set_env(dots)
}
