#'process_lfdr_results.R
#'
#' Sorts lfdrResults and calculates cumulative averages to return significant SNPs.
#'
#'@param orig_data data.frame with J rows, one for each composite null hypothesis being tested, holds data such as position, z-statistics, etc.
#'@param lfdrResults J*1 vector of all lfdr-values, the jth value corresponds to the jth row of orig_data.
#'@param fdr_threshold Scalar between 0 and 1 determining the percentage of false discoveries, default set to 0.1.
#'
#'@return A list with elements:
#'\item{all_snps_info}{A data.frame of all SNPs, sorted by the lfdr-value at which they are significant (cum_avg_lfdr column)}
#'\item{sig_snps_info}{A data.frame of the significant SNPs, sorted by the lfdr-value at which they are significant (cum_avg_lfdr column)}
#'\item{n_significant}{A scalar returning the number of significant SNPs.}
#'
#'@importFrom dplyr mutate filter arrange select
#'@importFrom rlang .data
#'
#'@export
#'
#'@examples
#'set.seed(0)
#'orig_data <- data.frame(
#'variant_id = sample(c("A","B","C","D"), 10^4, replace = TRUE),
#'Chr = sample(1:22, 10^4, replace = TRUE),
#'BP = sample(1e6:2e6, 10^4),
#'Z_Dataset1  = rnorm(10^4),
#'Z_Dataset2    = rnorm(10^4)
#')
#'lfdrResults <- runif(10^4)
#'process_lfdr_results(orig_data = orig_data, lfdrResults = lfdrResults, fdr_threshold = 0.05)


process_lfdr_results <- function(
    orig_data,
    lfdrResults,
    fdr_threshold = 0.1
) {

  #----------------------------------------------
  # Sort lfdr and compute cumulative averages
  #----------------------------------------------

  sorted_idx <- order(lfdrResults)

  cum_avg <- cumsum(
    lfdrResults[sorted_idx]
  ) / seq_along(sorted_idx)

  # map back to original SNP order
  cum_avg_original <- numeric(length(lfdrResults))

  cum_avg_original[sorted_idx] <- cum_avg

  #------------------------------------------------------
  # Attach lfdr results to merged data
  #------------------------------------------------------

  all_snps_info <- orig_data %>%
    mutate(
      lfdrResults = lfdrResults,
      cum_avg_lfdr = cum_avg_original
    ) %>%
    arrange(.data$cum_avg_lfdr)
  sig_snps_info <- all_snps_info %>% filter(.data$cum_avg_lfdr < fdr_threshold)

  #------------------------
  # Return results
  #------------------------

  return(list(
    all_snps_info = all_snps_info,
    sig_snps_info = sig_snps_info,
    n_significant = nrow(sig_snps_info)
  ))
}
