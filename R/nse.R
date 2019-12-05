#' Create a data mask rectangular data structures in Bioconductor
#' 
#' @description Tidy evaluation for rectangular 
#' data structures in Bioconductor.
#' @details This is the backend for non-standard evaluation 
#' in `plyranges` and is used to provide semantics for non-standard
#' evaluation used throughout the grammar. Generally,
#' you will not need to interact with this function directly,
#' but it can be useful if you're planning on extending
#' `plyranges` functionality. 
#' 
#' @seealso [rlang::new_data_mask()], [rlang::eval_tidy()]
#' @keywords internal
#' @rdname nse
#' @export
bc_data_mask <- function(data) UseMethod("bc_data_mask")

bc_data_mask.Vector <- function(data) {
  # extract the namespace of the class
  pkg_scope <- rlang::pkg_env(packageSlot(class(data)))
  
  top <- bc_fn_env()
  spec <- bc_fn_specials(data, top)
  mid <- bc_mcols_active(data, spec)
  bottom <- bc_vec_active(data, top, pkg_scope)
  
  mask <- rlang::new_data_mask(bottom, top = top)
  mask$.data <- rlang::as_data_pronoun(mask)
  mask
}

# extract generics and place them into an environment
bc_fn_env <- function(data) {
  top <- bioc_generics()
  rlang::new_environment(top)
}

# plyranges and dplyr special functions
bc_fn_specials <- function(data, env) {
  UseMethod("bc_fn_specials")
}

bc_fn_specials.Vector <- function(data, env) {
  rlang::env_bind_active(env, 
                         n = function() length(data))
  env
}

bc_fn_specials.List <- function(data, env) {
  special <- rlang::env(env)
  rlang::env_bind_active(special, 
                         n = function() lengths(data))
  special
}

bc_fn_specials.DFrame <- function(data, env) {
  special <- rlang::env(env)
  rlang::env_bind_active(special, 
                         n = function() nrow(data))
  special
}

bc_fn_specials.GroupedGenomicRanges <- function(data, env) {
  special <- rlang::env(env)
  rlang::env_bind_active(special, 
                         n = function() {
                           IRanges::runLength(dplyr::group_indices(data))
                           })
special
}


# 
bc_mcols_active <- function(data, env) {
  # enclose the mcols as middle
  mcols_names <- names(mcols(data))
  mcols_fn <- lapply(mcols_names,
                     function(nm) {
                       function() mcols(data)[[nm]]
                     })
  names(mcols_fn) <- mcols_names
  mid <- rlang::env(env)
  rlang::env_bind_active(mid, !!!mcols_fn)
  mid
}

bc_vec_active <- function(data, env, scope) {
  UseMethod("bc_vec_active")
}

bc_vec_active.Vector <- function(data, env, scope) {
  vec_names <- S4Vectors::parallelVectorNames(data)
  vec_fn <- lapply(vec_names,
                   function(nm) {
                     getter <- rlang::env_get(scope, nm)
                     function() getter(data)
                   })
  names(vec_fn) <- vec_names
  bottom <- rlang::env(env)
  rlang::env_bind_active(bottom, !!!vec_fn)
  bottom
}

bc_vec_active.List <- function(data, env, scope) {
  # bottom is the vector
  vec_names <- names(data)
  vec_fn <- lapply(vec_names,
                   function(nm) {
                     function() data[[nm]]
                   })
  names(vec_fn) <- vec_names
  bottom <- rlang::env(env)
  rlang::env_bind_active(bottom, !!!vec_fn)
  bottom
}



