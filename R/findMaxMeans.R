#' find_max_means.R
#'
#' Find maximum means for each dimension in null settings.

#' @param muInfo A list with 2^K elements, where each element is a matrix with K rows and Mb columns.
#'
#' @return A K*1 vector of the maximum means for each dimension under the null.
#'
#' @export
#' @examples
#' initMuList <- list(matrix(data=0, nrow=2, ncol=1), matrix(data=c(0, 3, 0, 6), nrow=2),
#' matrix(data=c(3, 0, 6, 0), nrow=2), matrix(data=c(8, 8), nrow=2))
#' find_max_means(initMuList)
find_max_means <- function(muInfo) {

  # iterate, skip the first (0) and last (alternative)
  listLength <- length(muInfo)
  K <- nrow(muInfo[[1]])
  # just keep finding the max
  maxMeans <- rep(0, K) 
  for (element_it in 2:(listLength - 1)) {
    tempMat <- cbind(muInfo[[element_it]], maxMeans)
    maxMeans <- apply(tempMat, 1, max)
  }
  # return K*1 vector
  return(maxMeans)
}





