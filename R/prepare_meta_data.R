prepare_meta_data <- function(
    datasets,
    dataset_names = NULL
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
    # Calculate Beta, SE, and Z-scores
    # -------------------------------
    beta_out <- paste0("Beta_", dataset_names[k])
    se_out   <- paste0("SE_", dataset_names[k])
    z_out <- paste0("Z_", dataset_names[k])
    p_out <- paste0("P_", dataset_names[k])
    
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
        !!beta_out := .data[[beta_col[1]]],
        !!se_out   := .data[[se_col[1]]],
        !!z_out := .data[[beta_col[1]]] / .data[[se_col[1]]],
        !!p_out := .data[[pval_col[1]]]
      ) %>%
      dplyr::filter(
        !is.na(.data[[z_out]]),
        !is.na(.data[[p_out]]),
        !is.na(.data[[beta_out]]) & !is.na(.data[[se_out]])
      ) %>%
      dplyr::select(dplyr::all_of(c("Chr", "BP", "A1", "A2", "rsid", beta_out, se_out, z_out, p_out))
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
        dplyr::mutate(across(all_of(c(beta_out, z_out)), \(x) if_else(.data$flip == 1, -x, x))) %>%
        dplyr::select(-dplyr::all_of(c("A1", "A2", "match", "flip")))
    }
    
  }
  
  # -----------------------------------
  # 
  # -----------------------------------
  betaDat <-
    merged_dat %>%
    select(starts_with("Beta_")) %>%
    as.matrix()
  
  seDat <-
    merged_dat %>%
    select(starts_with("SE_")) %>%
    as.matrix()
  
  
  return(list(
    extra_dat = merged_dat,
    betaDat = betaDat,
    seDat = seDat,
    n_snps = nrow(merged_dat)
  ))
}
