#' @useDynLib sdbmsABC
#' @importFrom Rcpp sourceCpp
NULL

#'@rdname Splitting_JRNMM_Cpp
#'@title Splitting_JRNMM_Cpp
#'@description Simulate a 6-dimensional trajectory from the JR-NMM
#' using Strang splitting
#'@return the trajectory of X=(X1,...,X6)
#'@export
Splitting_JRNMM_Cpp <- function(h, startv, grid, dm, cm, mu, C,A,B,a,b,v0,r,vmax){
  return(Splitting_JRNMM_Cpp_(h, startv, grid, dm, cm, mu, C,A,B,a,b,v0,r,vmax))
}

#'@rdname Splitting_JRNMM_output_Cpp
#'@title Splitting_JRNMM_output_Cpp
#'@description Simulate the output trajectory from the JR-NMM
#' using Strang splitting
#'@return the trajectory Y=X2-X3
#'@export
Splitting_JRNMM_output_Cpp <- function(h, startv, grid, dm, cm, mu, C,A,B,a,b,v0,r,vmax){
  return(Splitting_JRNMM_output_Cpp_(h, startv, grid, dm, cm, mu, C,A,B,a,b,v0,r,vmax))
}

#'@rdname Splitting_JRNMM_gen_Cpp
#'@title Splitting_JRNMM_gen_Cpp
#'@description Simulate a 6-dimensional trajectory from the JR-NMM
#' using Strang splitting when the matrix*vector multiplication is provided in a general form
#'@return the trajectory X=(X1,...,X6)
#'@export
Splitting_JRNMM_gen_Cpp <- function(h, startv, grid, dm, cm, mu, C,A,B,a,b,v0,r,vmax){
  return(Splitting_JRNMM_gen_Cpp_(h, startv, grid, dm, cm, mu, C,A,B,a,b,v0,r,vmax))
}
