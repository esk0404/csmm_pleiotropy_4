library(readr)
library(dplyr)
library(tidyr)
library(rtracklayer)
library(GenomicRanges)
library(data.table)
library(ggplot2)

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
    
    if (length(beta_col) == 0) {
      stop(paste0(dataset_names[k], "missing beta column"))
    }
    
    if (length(se_col) == 0) {
      stop(paste0(dataset_names[k], "missing standard error column"))
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
    # Calculate Beta and SE
    # -------------------------------
    beta_out <- paste0("Beta_", dataset_names[k])
    se_out   <- paste0("SE_", dataset_names[k])
    
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
        !!se_out   := .data[[se_col[1]]]
      ) %>%
      dplyr::filter(
        !is.na(.data[[beta_out]]) & !is.na(.data[[se_out]])
      ) %>%
      dplyr::select(dplyr::all_of(c("Chr", "BP", "A1", "A2", "rsid", beta_out, se_out))
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
      merged_dat <- merge(merged_dat, df, by=c("Chr", "BP")) %>%
        mutate(match = ifelse(.data$EA == .data$A1 & .data$OA == .data$A2, 1, 0)) %>%
        mutate(flip = ifelse(.data$EA == .data$A2 & .data$OA == .data$A1, 1, 0)) %>%
        dplyr::filter(.data$match + .data$flip > 0) %>%
        dplyr::mutate(across(all_of(beta_out), \(x) if_else(.data$flip == 1, -x, x))) %>%
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


# ----------------------------
# 0. Load data
# ----------------------------
setwd("/Users/ekim14/Downloads")

prostate <- read_tsv("GCST006085.txt",
                     col_select = c(Chr, position, Pvalue, rs_id, 
                                    SNP, Allele2, Allele1,
                                    Effect, StdErr),
                     progress = FALSE)

breast <- read_tsv("GCST004988.txt",
                   col_select = c(chr, var_name, phase3_1kg_id,
                                  position_b37, bcac_gwas_all_beta, 
                                  bcac_gwas_all_se,
                                  bcac_gwas_all_P1df,
                                  a0, a1),
                   progress = FALSE)

lung <- read_tsv("GCST004744_buildGRCh37.tsv",
                 col_select = c(variant_id, chromosome,
                                base_pair_location,
                                beta,
                                standard_error,
                                effect_allele,
                                other_allele,
                                p_value),
                 progress = FALSE)

ovarian <- read_tsv("GCST90455658.tsv.gz",
                    col_select = c(variant_id, chromosome, base_pair_location, 
                                   other_allele, effect_allele,
                                   beta, standard_error,
                                   p_value),
                    progress = FALSE)


prostate_new <- prostate %>%
  transmute(rsid = rs_id,
            chr = Chr,
            position = position,
            effect_allele = toupper(as.character(Allele1)),
            other_allele = toupper(as.character(Allele2)),
            pvalue = as.numeric(Pvalue),
            beta = Effect,
            se = StdErr
  ) %>%
  filter(!is.na(pvalue))

rm(prostate)

breast_new <- breast %>%
  transmute(rsid = phase3_1kg_id,
            chr = as.numeric(chr),
            position = position_b37,
            variant_id = var_name,
            pvalue = as.numeric(bcac_gwas_all_P1df),
            beta = as.numeric(bcac_gwas_all_beta),
            se = as.numeric(bcac_gwas_all_se),
            effect_allele = a0,
            other_allele = a1
  ) %>%
  filter(!is.na(pvalue))

rm(breast)
            
lung_new <- lung %>%
  transmute(rsid = variant_id,
            chr = chromosome,
            position = base_pair_location,
            effect_allele = effect_allele,
            other_allele = other_allele,
            beta = beta,
            se = standard_error,
            pvalue = p_value) %>%
  filter(!is.na(pvalue))

rm(lung)


chain_38to37 <- import.chain("hg38ToHg19.over.chain")

ovarian_liftover <- {
  keep0 <- !is.na(ovarian$base_pair_location) # get rid of NAs
  df <- ovarian[keep0, ]
  rm(ovarian); gc()
  
  df$chromosome <- as.character(df$chromosome)
  df$chromosome[df$chromosome == "23"] <- "X"
  
  gr <- GRanges(
    seqnames = paste0("chr", df$chromosome),
    ranges   = IRanges(
      start = df$base_pair_location,
      end   = df$base_pair_location
    )
  )
  
  genome(gr) <- "hg38"
  
  
  lifted <- liftOver(gr, chain_38to37)
  
  # keep only uniquely mapped variants
  keep1 <- elementNROWS(lifted) == 1
  lifted_gr <- unlist(lifted[keep1])
  
  df <- df[keep1, ]
  
  new_chr <- gsub("^chr", "", as.character(seqnames(lifted_gr)))
  
  df$chromosome <- as.numeric(new_chr)
  df$base_pair_location <- start(lifted_gr)
  
  df
}


rm(chain_38to37, keep0, df, gr, new_chr, keep1, lifted, lifted_gr)


datasets <- list(breast_new,lung_new,ovarian_liftover,prostate_new)
prepped_meta_data <- prepare_meta_data(datasets=datasets, dataset_names = c("Breast", "Lung", "Ovarian", "Prostate")) 







# Random effects meta-analysis
library(metafor)
library(dplyr)

run_random_meta <- function(prepped_meta_data) {
  
  betaDat <- prepped_meta_data$betaDat
  seDat   <- prepped_meta_data$seDat
  extra   <- prepped_meta_data$extra_dat
  
  # Remove SNPs with beta = 0 and SE = 0
  keep <- rowSums(betaDat == 0 & seDat == 0) == 0
  
  betaDat <- betaDat[keep, ]
  seDat   <- seDat[keep, ]
  extra   <- extra[keep, ]
  
  results <- lapply(seq_len(nrow(betaDat)), function(i) {
    
    fit <- rma(
      yi = betaDat[i, ],
      sei = seDat[i, ],
      method = "REML",
      control = list(maxiter = 100)
    )
    
    data.frame(
      beta_RE = as.numeric(fit$b),
      se_RE = fit$se,
      z_RE = fit$zval,
      p_RE = fit$pval,
      tau2 = fit$tau2,
      I2 = fit$I2,
      Q = fit$QE,
      Q_p = fit$QEp
    )
  })
  
  bind_cols(
    extra %>% select(Chr, BP, EA, OA),
    bind_rows(results)
  )
}


random_meta_results <- run_random_meta(
  prepped_meta_data
)


# Fixed effects meta analysis and manhattan plot

create_meta_plots <- function(
    plotData
) {
  
  # -------------------------------
  # Format data
  # -------------------------------
  plotData <- plotData %>%
    mutate(
      Chr = as.numeric(.data$Chr),
      BP = as.numeric(.data$BP),
      p_value = as.numeric(.data$p_val)
    ) %>%
    filter(
      !is.na(.data$Chr),
      !is.na(.data$BP),
      !is.na(.data$p_value)
    )
  
  
  # -------------------------------
  # Arrange data by chromosome
  # -------------------------------
  plotData <- plotData %>%
    arrange(.data$Chr)
  
  uniqueChrs <- sort(unique(plotData$Chr))
  
  
  # -------------------------------
  # Chromosome lengths
  # -------------------------------
  chrCounts <- plotData %>%
    group_by(.data$Chr) %>%
    summarise(
      chrLength = max(.data$BP, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(.data$Chr)
  
  
  # -------------------------------
  # Chromosome offsets
  # -------------------------------
  chrOffsets <- cumsum(chrCounts$chrLength)
  names(chrOffsets) <- chrCounts$Chr
  
  truePos <- numeric(nrow(plotData))
  
  for (chr_it in seq_along(uniqueChrs)) {
    
    tempChr <- uniqueChrs[chr_it]
    idx <- plotData$Chr == tempChr
    
    offsetVal <- if (chr_it == 1) {
      0
    } else {
      chrOffsets[chr_it - 1]
    }
    
    truePos[idx] <- plotData$BP[idx] + offsetVal
  }
  
  
  # -------------------------------
  # X-axis chromosome labels
  # -------------------------------
  xBreaks <- chrOffsets
  
  xBreaksLabs <- ifelse(
    uniqueChrs %% 2 == 0,
    "",
    as.character(uniqueChrs)
  )
  
  
  # -------------------------------
  # Add true genomic position
  # -------------------------------
  plotData <- plotData %>%
    mutate(
      truePos = truePos
    )
  
  
  # -------------------------------
  # Manhattan plot
  # -------------------------------
  
  y_max <- max(
    -log10(plotData$p_value),
    na.rm = TRUE
  ) + 1
  
  manPlot <- ggplot(
    plotData,
    aes(
      x = .data$truePos,
      y = -log10(.data$p_value)
    )
  ) +
    
    geom_point() +
    
    xlab("Chromosome") +
    ylab("-log10(p-value)") +
    
    scale_x_continuous(
      name = "Chr",
      breaks = xBreaks,
      labels = xBreaksLabs
    ) +
    
    ylim(
      c(0, y_max)
    ) +
    
    theme_cowplot()
  
  
  # -------------------------------
  # Z-score plot
  # -------------------------------
  
  zPlot <- ggplot(
    plotData,
    aes(
      x = .data$Z_Dataset1,
      y = .data$Z_Dataset2
    )
  ) +
    
    geom_point(
      color = "red"
    ) +
    
    geom_vline(
      xintercept = 0,
      linetype = "dashed"
    ) +
    
    geom_hline(
      yintercept = 0,
      linetype = "dashed"
    ) +
    
    xlab("Z-score 1") +
    ylab("Z-score 2") +
    
    theme_cowplot()
  
  
  # -------------------------------
  # Return plots
  # -------------------------------
  
  return(
    list(
      manPlot = manPlot,
      zPlot = zPlot
    )
  )
}

# Inverse-variance weighted (IVW) Fixed effects meta-analysis
B <- prepped_meta_data$betaDat
S <- prepped_meta_data$seDat

V <- S^2 # variance
W_FE <- 1 / V

beta_FE <- rowSums(B * W_FE, na.rm = TRUE) / rowSums(W_FE, na.rm = TRUE)
se_FE <- sqrt(1 / rowSums(W_FE, na.rm = TRUE))
z_FE <- beta_FE / se_FE
p_FE <- 2 * pnorm(-abs(z_FE))

fe_meta_df <- prepped_meta_data$extra_dat %>%
  mutate(
    beta_FE = beta_FE,
    se_FE   = se_FE,
    z_FE    = z_FE,
    p_val    = p_FE,
    variant_id = paste(
      Chr, BP, EA, OA,
      sep = "_"
    )
  ) %>%
  select(variant_id, Chr, BP, EA, OA, beta_FE, se_FE, z_FE, p_val)

fe_sig <- fe_meta_df %>%
  filter(p_val < 5e-8) 

manData_fe_t2 <- manData_t2 %>%
  mutate(
    FE_sig = paste(Chr, BP, sep = "_") %in%
      paste(fe_sig$Chr, fe_sig$BP, sep = "_")
  )


figs <- create_meta_plots(plotData = fe_sig)
figs$manPlot


## Save plot
ggsave(
  filename = "meta_manplot.png",
  plot = figs$manPlot,
  width = 17,
  height = 7,
  dpi = 300
)




