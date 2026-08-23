#' calc_dens_ind.R
#'
#' Calculate J bivariate normal densities (both dimensions are independent) under fitted csmGmm.
#'
#' @param x 2*1 vector of means.
#' @param Zmat J*2 matrix of test statistics.
#'
#' @return A J*1 vector of densities for each row of Zmat.
#' @importFrom stats dnorm
#'
#' @export
#' @examples
#' x <- c(0, 0)
#' Zmat <- cbind(rnorm(10^5), rnorm(10^5))
#' calc_dens_ind_2d(x, Zmat)
#'
calc_dens_ind_2d <- function(x, Zmat) {
  # right now, only for 2 or 3 dimensions
  stats::dnorm(Zmat[, 1], mean=x[1], sd=1) * stats::dnorm(Zmat[, 2], mean=x[2], sd=1)
}


#' Calculate J trivariate normal densities (all dimensions are independent) under fitted csmGmm.
#'
#' @param x 3*1 vector of means.
#' @param Zmat J*3 matrix of test statistics.
#'
#' @return A J*1 vector of densities for each row of Zmat.
#' @importFrom stats dnorm
#'
#' @export
#' @examples
#' x <- c(0, 0)
#' Zmat <- cbind(rnorm(10^5), rnorm(10^5), rnorm(10^5))
#' calc_dens_ind_3d(x, Zmat)
#'
calc_dens_ind_3d <- function(x, Zmat) {
  # right now, only for 2 or 3 dimensions
  stats::dnorm(Zmat[, 1], mean=x[1], sd=1) * stats::dnorm(Zmat[, 2], mean=x[2], sd=1) * stats::dnorm(Zmat[, 3], mean=x[3], sd=1)
}



#' Calculate the density of K-dimensional multivariate normal (all dimensions are independent) under fitted acsGmm.
#'
#' @param x K*1 vector of means.
#' @param Zmat J*K matrix of test statistics.
#'
#' @return A J*1 vector of densities for each row of Zmat.
#' @importFrom stats dnorm
#'
#' @export
#' @examples
#' x <- c(0, 0)
#' Zmat <- cbind(rnorm(10^5), rnorm(10^5), rnorm(10^5), rnorm(10^5))
#' calc_dens_ind_multiple(x, Zmat)
#'
calc_dens_ind_multiple <- function(x, Zmat) {
  K <- ncol(Zmat)
  tempSum <- rep(0, nrow(Zmat))
  for (k_it in 1:K) {
    tempSum <- tempSum + stats::dnorm(Zmat[, k_it], mean=x[k_it], sd=1, log=TRUE)
  }
  exp(tempSum)
}
