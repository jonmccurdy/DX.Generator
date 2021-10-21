# this can be used for other parameters.
#' set.dx
#'
#' @param K parameter for dx generator
#' @param S parameter for dx generator
#'
#' @return Creates an initialized table of size k
#' @export
#'
#' @examples set.dx(K=47, S=1)
set.dx <- function(K=47, S=1){
  if (RNGkind()[1] != "user-supplied"){RNGkind("user-supplied")}
  if (K>623){
    message('The value of K must be less than 624. By default, K=47 chosen.')
    K=47
  }
  set.seed(.Random.seed[4])
}
