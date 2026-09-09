# ============================================================
# Run ASSET meta-analysis
# ============================================================

# Load packages
library(dplyr)
library(rtracklayer)
library(GenomicRanges)
library(data.table)
library(ASSET)
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
# 6. Define sample sizes
# ============================================================

ncase <- c(
  Breast = 76192,
  Lung = 11273,
  Ovarian = 15361,
  Prostate = 79148
)

ncntl <- c(
  Breast = 63082,
  Lung = 55483,
  Ovarian = 104887,
  Prostate = 61106
)


# ============================================================
# 7. Define correlation matrix
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


# ============================================================
# 8. Run ASSET
# ============================================================

run_ASSET <- function(
    prepped_meta_data,
    ncase,
    ncntl,
    cor_mat
) {

  beta <- prepped_meta_data$betaDat
  se <- prepped_meta_data$seDat
  extra <- prepped_meta_data$extra_dat

  # ASSET requires at least 3 valid cancer-specific estimates
  keep <- rowSums(
    is.finite(beta) &
      is.finite(se) &
      se > 0
  ) >= 3

  beta <- beta[keep, , drop = FALSE]
  se <- se[keep, , drop = FALSE]
  extra <- extra[keep, , drop = FALSE]


  # ----------------------------------------------------------
  # One-sided subset test
  # ----------------------------------------------------------

  asset_one <- h.traits(
    snp.vars = extra$rsid,
    traits.lab = c(
      "Breast",
      "Lung",
      "Ovarian",
      "Prostate"
    ),
    beta.hat = beta,
    sigma.hat = se,
    ncase = ncase,
    ncntl = ncntl,
    cor = cor_mat,
    search = 1,
    side = 2,
    meta = TRUE,
    meth.pval = "DLM"
  )


  # ----------------------------------------------------------
  # Two-sided subset test
  # ----------------------------------------------------------

  asset_two <- h.traits(
    snp.vars = extra$rsid,
    traits.lab = c(
      "Breast",
      "Lung",
      "Ovarian",
      "Prostate"
    ),
    beta.hat = beta,
    sigma.hat = se,
    ncase = ncase,
    ncntl = ncntl,
    cor = cor_mat,
    search = 2,
    meta = TRUE,
    meth.pval = "DLM"
  )


  # ----------------------------------------------------------
  # Extract ASSET results
  # ----------------------------------------------------------

  one_sided_p <-
    asset_one$Subset.1sided$pval

  one_fe_p <-
    asset_one$Meta$pval

  two_sided_p <-
    asset_two$Subset.2sided$pval

  two_sided_pos_p <-
    asset_two$Subset.2sided$pval.1

  two_sided_neg_p <-
    asset_two$Subset.2sided$pval.2


  # ----------------------------------------------------------
  # Combine results
  # ----------------------------------------------------------

  asset_results <- extra %>%
    mutate(
      asset1_sided_pval = one_sided_p,
      asset2_sided_pval = two_sided_p,
      asset2_sided_pval_pos = two_sided_pos_p,
      asset2_sided_pval_neg = two_sided_neg_p,
      asset_fe_pval = one_fe_p
    ) %>%
    select(
      "Chr",
      "BP",
      "EA",
      "OA",
      "rsid",
      asset1_sided_pval,
      asset2_sided_pval,
      asset2_sided_pval_pos,
      asset2_sided_pval_neg,
      asset_fe_pval
    )


  return(asset_results)
}


asset_results <- run_ASSET(
  prepped_meta_data = prepped_meta_data,
  ncase = ncase,
  ncntl = ncntl,
  cor_mat = cor_mat
)


# ============================================================
# 9. Filter and sort significant results
# ============================================================

asset1_sig_results <- asset_results %>%
  filter(asset1_sided_pval < 1.25e-8) %>%
  arrange(asset1_sided_pval)


asset2_sig_results <- asset_results %>%
  filter(asset2_sided_pval < 1.25e-8) %>%
  arrange(asset2_sided_pval)


# ============================================================
# 10. Save results
# ============================================================

write.csv(
  asset_results,
  file.path(
    output_dir,
    "asset_results.csv"
  ),
  row.names = FALSE
)


write.csv(
  asset1_sig_results,
  file.path(
    output_dir,
    "asset1_sided_sig.csv"
  ),
  row.names = FALSE
)


write.csv(
  asset2_sig_results,
  file.path(
    output_dir,
    "asset2_sided_sig.csv"
  ),
  row.names = FALSE
)
