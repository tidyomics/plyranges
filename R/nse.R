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
bc_data_mask <- function(data) {
  # extract the namespace of the class
  pkg_scope <- rlang::pkg_env(packageSlot(class(data)))
  
  top <- bc_fn_env(data)
  mid <- bc_vec_active(mcols(data), top)
  bottom <- bc_vec_active(data, mid, pkg_scope)
  
  mask <- rlang::new_data_mask(bottom, top = top)
  mask$.data <- rlang::as_data_pronoun(mask)
  mask
}

# extract generics and place them into an environment
bc_fn_env <- function(data) {
  top <- bioc_generics()
  top <- rlang::new_environment(top)
  bc_fn_specials(data, top)
}

# plyranges and dplyr special functions
bc_fn_specials <- function(data, env) {
  UseMethod("bc_fn_specials")
}

bc_fn_specials.default <- function(data, env) {
  stopifnot(is(data, "Vector"))
  rlang::env_bind_active(env,
                         n = ~ function() length(data))
  env
}

bc_fn_specials.DFrame <- function(data, env) {
  rlang::env_bind_active(env, 
                         n = ~ function() nrow(data))
  env
}

bc_fn_specials.GroupedGenomicRanges <- function(data, env) {
  rlang::env_bind_active(env, 
                         n = ~ function() {
                           IRanges::runLength(dplyr::group_indices(data))
                         })
  env
}

bc_fn_specials.GroupedIntegerRanges <- bc_fn_specials.GroupedGenomicRanges

bc_vec_active <- function(data, env, ...) {
  UseMethod("bc_vec_active")
}

bc_vec_active.DFrame <- function(data, env) {
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


# eval_mask dispatch? 
bc_eval_tidy <- function(dots, data, mask = bc_data_mask(data)) {
  container <- DataFrame()
}


overscope_eval_update <- function(overscope, dots, bind_envir = TRUE) {
  update <- vector("list", length(dots))
  names(update) <- names(dots)
  for (i in seq_along(update)) {
    quo <- dots[[i]]
    update[[i]] <- eval_tidy(quo, data = overscope)
    # sometimes we want to compute on previously constructed columns
    # we can do this by binding the evaluated expression to
    # the overscope environment
    if (bind_envir) {
      new_col <- names(dots)[[i]]
      rlang::env_bind(overscope, !!new_col := update[[i]])
    }
  }
  return(update)
}





