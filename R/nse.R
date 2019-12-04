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
bc_data_mask <- function(data) {
  # extract the namespace of the class
  pkg_scope <- rlang::pkg_env(packageSlot(class(data)))
  
  # place the generics at the top of the mask
  top <- bioc_generics()
  top <- rlang::new_environment(top)
  # enclose the mcols as middle
  mcols_names <- names(mcols(data))
  mcols_fn <- lapply(mcols_names,
                   function(nm) {
                     function() mcols(data)[[nm]]
                   })
  names(mcols_fn) <- mcols_names
  
  mid <- rlang::env(top)
  rlang::env_bind_active(mid, !!!mcols_fn)
  # bottom is the vector
  vec_names <- parallelVectorNames(data)
  vec_fn <- lapply(vec_names,
                   function(nm) {
                     getter <- rlang::env_get(pkg_scope, nm)
                     function() getter(data)
                   })
  names(vec_fn) <- vec_names
  bottom <- rlang::env(mid)
  rlang::env_bind_active(bottom, !!!vec_fn)
  
  mask <- rlang::new_data_mask(bottom, top = top)
  mask$.data <- rlang::as_data_pronoun(mask)
  mask
}

