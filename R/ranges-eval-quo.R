# A half baked algortihm for statically parsing a call from a quosure
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


