#' prepare_csmgmm_data.R
#'
#' Prepares J*K matrix of summary statistics for GWAS replication studies where J is the number of SNPs and K is the number of cohorts.
#'
#' @param datasets list holding K datasets to use for replication analysis
#' @param dataset_names optional vector of names for each dataset; otherwise function uses Dataset1,.., DatasetK as default.
#' @param z_cap optional number to cap the most extreme z-scores, we recommend setting at 8.1 as R is not very accurate for very large z-scores
#'
#' @return A list with the elements:
#' \item{extra_dat}{A J*K matrix with J matching SNPs across K datasets, holds relevant SNP information.}
#' \item{clean_dat}{J*K matrix of test statistics, holds J Z-scores for each dataset.}
#' \item{n_snps}{A scalar returning the number of matching SNPs across all K datasets.}
#'
#' @importFrom dplyr %>% mutate filter select across if_else all_of starts_with
#' @importFrom rlang .data :=
#'
#' @export
#' @examples
#' library(dplyr)
#' set.seed(0)
#' J <- 1000
#' alleles <- c("A", "C", "T", "G")
#' dat1 <- data.frame(zstat = rnorm(J)) %>% mutate(se = runif(n=J, min=0.0005, max=0.0015)) %>%
#' mutate(beta = zstat * se, chr=1, pos=1:J, a1=sample(x=alleles, J, replace=TRUE),
#' a2=sample(x=alleles, J, replace=TRUE))
#' dat2 <- data.frame(zstat = rnorm(J)) %>% mutate(se = runif(n=J, min=0.0005, max=0.0015)) %>%
#' mutate(beta = zstat * se, chr=1, pos=1:J, a1=sample(x=alleles, J, replace=TRUE),
#' a2=sample(x=alleles, J, replace=TRUE))
#' dat3 <- data.frame(zstat = rnorm(J)) %>% mutate(se = runif(n=J, min=0.0005, max=0.0015)) %>%
#' mutate(beta = zstat * se, chr=1, pos=1:J, a1=sample(x=alleles, J, replace=TRUE),
#' a2=sample(x=alleles, J, replace=TRUE))
#' datasets = list(dat1, dat2, dat3)
#' prepare_csmgmm_data(datasets = datasets)
#'
prepare_csmgmm_data <- function(
    datasets,
    dataset_names = NULL,
    z_cap = 8.1
) {
  
  # -----------------------------------
  # Dataset names
  # -----------------------------------
  if (is.null(dataset_names)) {
    dataset_names <- paste0("Dataset", seq_along(datasets))
  }
  
  if (length(dataset_names) != length(datasets)) {
    stop("dataset_names must match datasets")
  }
  
  merged_dat <- c()
  for (k in 1:length(datasets)) {
    
    # -----------------------------------------------
    # Check required columns for Z-score construction
    # -----------------------------------------------
    beta_col <- grep("(beta|^effect$|effect(?!_?allele))", names(datasets[[k]]), value = TRUE, ignore.case = TRUE, perl = TRUE)
    se_col <- grep("(std[_]?err|standard[_]?error|se$)", names(datasets[[k]]), value = TRUE, ignore.case = TRUE)
    pval_col <- grep("(pval|p_value|pvalue|p[_]?value)", names(datasets[[k]]), value = TRUE, ignore.case = TRUE)
    
    
    if (length(beta_col) == 0) {
      stop(paste0(dataset_names[k], "missing beta column"))
    }
    
    if (length(se_col) == 0) {
      stop(paste0(dataset_names[k], "missing standard error column"))
    }
    
    if (length(pval_col) == 0) {
      stop(paste0(dataset_names[k], " missing p-value column"))
    }
    
    
    # -------------------------------
    # Standardize variant information
    # -------------------------------
    chr_col <- grep("(^chr|chromosome)", names(datasets[[k]]), value = TRUE, ignore.case = TRUE)[1]
    position_col <- grep("(pos|loc|position|bp)", names(datasets[[k]]), value = TRUE, ignore.case = TRUE)[1]
    ea_col <- grep("(^a.*1|effect_allele|\\bea\\b)", names(datasets[[k]]), value = TRUE, ignore.case = TRUE)[1]
    oa_col <- grep("(^a.*2|other_allele|\\boa\\b)", names(datasets[[k]]), value = TRUE, ignore.case = TRUE)[1]
    rsid_col <- grep("(rsid|rs_id)", names(datasets[[k]]), value = TRUE, ignore.case = TRUE)[1]
    required_id_cols <- c(chr_col, position_col, oa_col, ea_col)
    if (any(is.na(required_id_cols))) {
      stop(paste0(dataset_names[k], " cannot create variant ID (missing columns)"))
    }
    
    # -------------------------------
    # Calculate Z-scores and add variant ID
    # -------------------------------
    z_col_name <- paste0("Z_", dataset_names[k])
    p_col_name <- paste0("P_", dataset_names[k])
    
    df <- datasets[[k]] %>%
      mutate(
        Chr = .data[[chr_col]],
        A1 = .data[[ea_col]],
        A2 = .data[[oa_col]],
        BP = .data[[position_col]],
        rsid = if (!is.na(rsid_col)) {
          .data[[rsid_col]]
        } else {
          NA_character_
        },
        !!z_col_name := .data[[beta_col[1]]] / .data[[se_col[1]]],
        !!p_col_name := .data[[pval_col[1]]]
      ) %>%
      dplyr::filter(
        !is.na(.data[[z_col_name]]),
        !is.na(.data[[p_col_name]])
      ) %>%
      dplyr::select(dplyr::all_of(c("Chr", "BP", "A1", "A2", "rsid", z_col_name, p_col_name))
      ) %>%
      dplyr::filter(
        nchar(paste0(.data$A1, .data$A2)) == 2
      )
    
    
    # -----------------------------------------------
    # Flip for allele if needed
    # -----------------------------------------------
    if (k == 1) {
      merged_dat <- df
      colnames(merged_dat)[3:4] <- c("EA", "OA")
    } else {
      df <- df %>% select(-rsid)
      
      merged_dat <- merge(merged_dat, df, by=c("Chr", "BP")) %>%
        mutate(match = ifelse(.data$EA == .data$A1 & .data$OA == .data$A2, 1, 0)) %>%
        mutate(flip = ifelse(.data$EA == .data$A2 & .data$OA == .data$A1, 1, 0)) %>%
        dplyr::filter(.data$match + .data$flip > 0) %>%
        dplyr::mutate(across(all_of(z_col_name), \(x) if_else(.data$flip == 1, -x, x))) %>%
        dplyr::select(-dplyr::all_of(c("A1", "A2", "match", "flip")))
    }
    
  }
  
  # -----------------------------------
  # Extract Z-score matrix
  # -----------------------------------
  testDat <- merged_dat %>%
    select(starts_with("Z_")) %>%
    as.matrix()
  
  # -----------------------------------
  # Cap Z-scores
  # -----------------------------------
  testDat[testDat > z_cap] <- z_cap
  testDat[testDat < -z_cap] <- -z_cap
  
  return(list(
    extra_dat = merged_dat,
    clean_dat = testDat,
    n_snps = nrow(merged_dat)
  ))
}
