#' @title dx_init
#' @name dx_init
#'
#' @description 
#' \code{dx_init()} allows the user to designate the DX-class of Pseudo-Random Number Generators 
#' as the user-supplied RNG in R. This makes it so where \code{runif()} and any other function 
#' relying on \code{runif()} will use a DX-k-s RNG to generate the random variate. 
#' 
#' @details 
#' Within this function there are two arguments which can be passed through the function.
#' \eqn{K} determines the order of the linear occurrence, and \eqn{S} determines the number of non-zero 
#' coefficients in the MRG. 
#' 
#' Speed is the fastest for the default values of \eqn{S}, that is \eqn{S = 1}.  
#' The argument \eqn{K} has no affect on the speed of the generator. The argument \eqn{S} essentially 
#' controls how many additions are performed to produce the next random number: the fewer, the faster.
#' 
#' Although DX-\eqn{k}-\eqn{s} generators have excellent empirical properties and can pass the most 
#' stringent battery of tests for uniformity, those set with \eqn{s = 4} will appear marginally more 
#' randomly and uniformly distributed. 
#' 
#' Often times in simulations, a consecutive run of random numbers need to appear independent with 
#' small serial correlations between any successive subset of numbers.  If \eqn{j} consecutive 
#' numbers need to appear independent, then choose argument \eqn{K} such that \eqn{K > j}.
#' 
#' @param K the order of the linear recurrence. \eqn{K} must be an integer between 5 and 50,000
#' @param S number of non-zero coefficients in the MRG. \eqn{S} must be an either 1,2,3 or 4
#'
#' @examples
#' \dontrun{dx_init(K=53, S=2)}
#' 
#' @references 
#' L.-Y. Deng and H. Xu. A system of high-dimensional, efficient, long-cycle and portable uniform random number generators. ACM Transactions on Modeling and Computer Simulation (TOMACS), 13 (4):299–309, 2003.

dx_init <- function(K=1597, S=1) {
  invisible(capture.output(setRNG::setRNG(kind="user-supplied", seed=c(K,S,.Random.seed[4])),type = "output"))
  invisible(runif(1))
  set.seed(.Random.seed[4])
}