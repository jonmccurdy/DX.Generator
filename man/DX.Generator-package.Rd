\name{DX.Generator-package}
\alias{DX.Generator}
\docType{package}
\title{
  DX.Generator implements the DX-class of RNG as a user-supplied generator in R, accessible with \code{runif()}
}
\description{
This package implements the DX-class of Psuedo-Random Number Generators which are very efficient and have nice empirical properties. When it comes to Random Number Generators, R currently provides seven different generators which are built-in; the default being the MT-19937 generator which has not passed all empirical tests. Users are able to select the DX-class of generators, which has passed all empirical test, as the user-supplied generator as well as choose the order of the linear recurrence as well as number of non-zero coefficients in the MRG by using the function.
}
\details{
  Users can set the DX-class of RNGs as the default generator in R by using the dx_init() function.
  Additionally, \code{dx_init(K=47, S=1)} will allow users to specify the K and S values for the generator.
}
\author{
Jonathan McCurdy, <jrmccrdy@memphis.edu>

Author: Lih-Yuan Deng
}
\references{
   L.-Y. Deng and H. Xu. A system of high-dimensional, efficient, long-cycle and portable uniform random number generators. ACM Transactions on Modeling and Computer Simulation (TOMACS), 13 (4):299–309, 2003.}
\keyword{ package }

\examples{
  \dontrun{
     dx_init()
     dx_init(K=53, S=2)
  }
}

