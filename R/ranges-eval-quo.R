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


# A half baked algortihm for dynamicall parsing a call from a quosure
# and mapping it over list objects
quo_set_generic <- function(quo, generics = bioc_generics()) {
  if (!rlang::quo_is_call(quo)) {
    return(quo)
  } else {
      fn_name <- rlang::call_name(quo)
      args_list <- rlang::call_args(quo)
      quo_env <- rlang::quo_get_env(quo)
      # idea here is to check for whether the function being called in the
      # quosure has has a method (with appropriate signature) 
      # in the quosoure environment
      # if it does we just return the quo
      # if it doesn't we add the generics from the bioc to quosoure environment
      # and repeat the process 
      # if there is still no method and if the signature is a list we modify
      # the call to pass to a Map 
  }
}