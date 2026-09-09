# ============================================================
# Run random-effects meta-analysis
# ============================================================

# Load packages
library(dplyr)
library(rtracklayer)
library(GenomicRanges)
library(data.table)
library(remaCor)
library(here)


# Load project-specific functions
source(here("R", "functions", "prepare_meta_data.R"))


# ============================================================
# 1. Set input/output directories
# ============================================================

data_dir <- here("data", "raw")
reference_dir <- here("data", "reference")
output_dir <- here("results", "meta_analysis")

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 2. Read GWAS summary statistics
# ============================================================

prostate_read <- fread(
  file.path(data_dir, "GCST006085.txt")
)

prostate <- prostate_read %>%
  select(
    Chr, position, Pvalue, rs_id,
    SNP, Allele2, Allele1,
    Effect, StdErr
  )


breast_read <- fread(
  file.path(data_dir, "GCST004988.txt")
)

breast <- breast_read %>%
  select(
    chr, var_name, phase3_1kg_id,
    position_b37, bcac_gwas_all_beta,
    bcac_gwas_all_se,
    bcac_gwas_all_P1df,
    a0, a1
  )


lung_read <- fread(
  file.path(data_dir, "GCST004744_buildGRCh37.tsv")
)

lung <- lung_read %>%
  select(
    variant_id, chromosome,
    base_pair_location,
    beta,
    standard_error,
    effect_allele,
    other_allele,
    p_value
  )


ovarian_read <- fread(
  file.path(data_dir, "GCST90455658.tsv.gz")
)

ovarian <- ovarian_read %>%
  select(
    variant_id,
    chromosome,
    base_pair_location,
    other_allele,
    effect_allele,
    beta,
    standard_error,
    p_value
  )


# ============================================================
# 3. Standardize datasets
# ============================================================

prostate_new <- prostate %>%
  transmute(
    rsid = rs_id,
    chr = Chr,
    position = position,
    effect_allele = toupper(as.character(Allele1)),
    other_allele = toupper(as.character(Allele2)),
    pvalue = as.numeric(Pvalue),
    beta = Effect,
    se = StdErr
  ) %>%
  filter(!is.na(pvalue))


breast_new <- breast %>%
  transmute(
    rsid = phase3_1kg_id,
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


lung_new <- lung %>%
  transmute(
    rsid = variant_id,
    chr = chromosome,
    position = base_pair_location,
    effect_allele = effect_allele,
    other_allele = other_allele,
    beta = beta,
    se = standard_error,
    pvalue = p_value
  ) %>%
  filter(!is.na(pvalue))


# ============================================================
# 4. Lift over ovarian cancer data: GRCh38 -> GRCh37
# ============================================================

chain_38to37 <- import.chain(
  file.path(
    reference_dir,
    "hg38ToHg19.over.chain"
  )
)


ovarian_liftover <- {

  keep0 <- !is.na(ovarian$base_pair_location)

  df <- ovarian[keep0, ]

  df$chromosome <- as.character(df$chromosome)

  df$chromosome[
    df$chromosome == "23"
  ] <- "X"


  gr <- GRanges(
    seqnames = paste0(
      "chr",
      df$chromosome
    ),
    ranges = IRanges(
      start = df$base_pair_location,
      end = df$base_pair_location
    )
  )

  genome(gr) <- "hg38"


  lifted <- liftOver(
    gr,
    chain_38to37
  )


  # Keep only uniquely mapped variants
  keep1 <- elementNROWS(lifted) == 1

  lifted_gr <- unlist(
    lifted[keep1]
  )

  df <- df[keep1, ]


  new_chr <- gsub(
    "^chr",
    "",
    as.character(
      seqnames(lifted_gr)
    )
  )

  new_chr[new_chr == "X"] <- "23"


  df$chromosome <- as.numeric(new_chr)

  df$base_pair_location <-
    start(lifted_gr)

  df
}


# ============================================================
# 5. Prepare meta-analysis data
# ============================================================

datasets <- list(
  breast_new,
  lung_new,
  ovarian_liftover,
  prostate_new
)

prepped_meta_data <- prepare_meta_data(
  datasets = datasets,
  dataset_names = c(
    "Breast",
    "Lung",
    "Ovarian",
    "Prostate"
  )
)


# ============================================================
# 6. Define correlation matrix
# ============================================================


cor_mat <- matrix(
  c(
    1.000,  0.031,  0.051,  0.004,
    0.031,  1.000,  0.013, -0.010,
    0.051,  0.013,  1.000, -0.005,
    0.004, -0.010, -0.005,  1.000
  ),
  nrow = 4,
  byrow = TRUE
)

rownames(cor_mat) <- colnames(cor_mat) <- c(
  "Breast",
  "Lung",
  "Ovarian",
  "Prostate"
)


run_RE2C <- function(prepped_meta_data, cor_mat) {
  
  beta <- prepped_meta_data$betaDat
  se <- prepped_meta_data$seDat
  extra <- prepped_meta_data$extra_dat
  
  keep <- rowSums(
    !is.finite(beta) |
      !is.finite(se) |
      se <= 0
  ) == 0
  
  beta <- beta[keep, , drop = FALSE]
  se <- se[keep, , drop = FALSE]
  extra <- extra[keep, , drop = FALSE]
  
  n <- nrow(beta)
  
  results <- vector("list", n)
  
  for (i in seq_len(n)) {
    
    results[[i]] <- RE2C(
      beta = beta[i, ],
      stders = se[i, ],
      cor = cor_mat
    )
  }
  
  re_results <- do.call(
    rbind,
    results
  )
  
  cbind(
    extra[, c("rsid", "Chr", "BP", "EA", "OA",
      "Z_Breast",
      "Z_Lung",
      "Z_Ovarian",
      "Z_Prostate",
      "P_Breast",
      "P_Lung",
      "P_Ovarian",
      "P_Prostate")],
    re_results
  )
}

re_results <- run_RE2C(
  prepped_meta_data,
  cor_mat
)



# ============================================================
# 7. Filter and sort significant results
# ============================================================

final_results <- re_results %>%
  filter(RE2Cp < 1.25e-8) %>%
  arrange(RE2Cp)


# ============================================================
# 8. Save results
# ============================================================

outfile <- file.path(
  output_dir,
  "re_meta_results.csv"
)

write.csv(
  final_results,
  outfile,
  row.names = FALSE
)
