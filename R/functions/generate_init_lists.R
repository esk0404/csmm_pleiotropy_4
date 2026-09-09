#' generate_init_lists.R
#'
#' Generates initMuList and initPiList parameters for symm_fit_ind_EM() function.
#'
#' @param K Scalar, number of datasets in the analysis.
#'
#' @return A list with the elements:
#' \item{initMuList}{List of 2^K elements where each element is a matrix with K rows and number of columns equal to the number of possible mean vectors for that binary representation}
#' \item{initPiList}{List of 2^K elements where each element is a vector with number of elements equal to the number of possible mean vectors for that binary representation.}
#' @export
#' @examples
#' K <- 2
#' generate_init_lists(K)



generate_init_lists <- function(K) {
  if (K == 2) {
    initPiList <- list(c(0.82), c(0.02, 0.02), c(0.02, 0.02), c(0.1))
    initMuList <- list(matrix(data=0, nrow=2, ncol=1), matrix(data=c(0, 3, 0, 6), nrow=2),
                       matrix(data=c(3, 0, 6, 0), nrow=2), matrix(data=c(8, 8), nrow=2))
  } else if (K == 3) {
    initPiList <- list(c(0.82))

    for (i in 2:7) {
      initPiList[[i]] <- c(0.08 / 12, 0.08 / 12)
    }

    initPiList[[8]] <- c(0.1)

    initMuList <- list(
      matrix(data = 0, nrow = 3, ncol = 1)
    )

    for (i in 2:7) {
      initMuList[[i]] <- cbind(rep(2, 3),
                               rep(5, 3))
    }

    initMuList[[8]] <- matrix(
      data = c(8, 8, 8),
      nrow = 3
    )

  } else if (K >= 4) {

    initPiList <- list(c(0.82))

    for (i in 2:(2^K)) {
      initPiList[[i]] <- 0.18 / (2^K - 1)
    }

    initMuList <- list(
      matrix(data = rep(0, K),
             nrow = K,
             ncol = 1)
    )

    for (i in 2:(2^K)) {
      initMuList[[i]] <- matrix(
        data = rep(3, K),
        nrow = K,
        ncol = 1
      )
    }

  } else {
    stop("K must be >= 2")
  }

  return(list(
    initPiList = initPiList,
    initMuList = initMuList
  ))
}
