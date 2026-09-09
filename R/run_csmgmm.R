# ============================================================
# Run csmGmm K=4 pleiotropy analysis
# ============================================================

# Load packages
library(csmGmm)
library(here)
library(readr)
library(dplyr)
library(GenomicRanges)
library(rtracklayer)


# Load project-specific functions
source(here("R", "functions", "prepare_csmgmm_data.R"))
source(here("R", "functions", "generate_init_lists.R"))
source(here("R", "functions", "process_lfdr_results.R"))
source(here("R", "functions", "symm_fit_ind_EM_R1.R"))


# ============================================================
# 1. Set input/output directories
# ============================================================

data_dir <- here("data", "raw")
reference_dir <- here("data", "reference")
output_dir <- here("results", "csmGmm")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# ============================================================
# 2. Read GWAS summary statistics
# ============================================================

prostate <- read_tsv(
  file.path(data_dir, "GCST006085.txt"),
  col_select = c(
    Chr, position, Pvalue, rs_id,
    SNP, Allele2, Allele1,
    Effect, StdErr
  ),
  progress = FALSE
)

breast <- read_tsv(
  file.path(data_dir, "GCST004988.txt"),
  col_select = c(
    chr, var_name, phase3_1kg_id,
    position_b37, bcac_gwas_all_beta,
    bcac_gwas_all_se,
    bcac_gwas_all_P1df,
    a0, a1
  ),
  progress = FALSE
)

lung <- read_tsv(
  file.path(data_dir, "GCST004744_buildGRCh37.tsv"),
  col_select = c(
    variant_id, chromosome,
    base_pair_location,
    beta, standard_error,
    effect_allele, other_allele,
    p_value
  ),
  progress = FALSE
)

ovarian <- read_tsv(
  file.path(data_dir, "GCST90455658.tsv.gz"),
  col_select = c(
    variant_id, chromosome,
    base_pair_location,
    other_allele, effect_allele,
    beta, standard_error,
    p_value
  ),
  progress = FALSE
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
  file.path(reference_dir, "hg38ToHg19.over.chain")
)

ovarian_liftover <- {
  
  keep0 <- !is.na(ovarian$base_pair_location)
  df <- ovarian[keep0, ]
  
  df$chromosome <- as.character(df$chromosome)
  df$chromosome[df$chromosome == "23"] <- "X"
  
  gr <- GRanges(
    seqnames = paste0("chr", df$chromosome),
    ranges = IRanges(
      start = df$base_pair_location,
      end = df$base_pair_location
    )
  )
  
  genome(gr) <- "hg38"
  
  lifted <- liftOver(gr, chain_38to37)
  
  # Keep only uniquely mapped variants
  keep1 <- elementNROWS(lifted) == 1
  lifted_gr <- unlist(lifted[keep1])
  
  df <- df[keep1, ]
  
  new_chr <- gsub(
    "^chr", "",
    as.character(seqnames(lifted_gr))
  )
  
  new_chr[new_chr == "X"] <- "23"
  
  df$chromosome <- as.numeric(new_chr)
  df$base_pair_location <- start(lifted_gr)
  
  df
}


# ============================================================
# 5. Prepare data for csmGmm
# ============================================================

datasets <- list(
  breast_new,
  lung_new,
  ovarian_liftover,
  prostate_new
)

prepped_data <- prepare_csmgmm_data(
  datasets = datasets,
  dataset_names = c(
    "Breast",
    "Lung",
    "Ovarian",
    "Prostate"
  )
)


# ============================================================
# 6. Initialize csmGmm parameters
# ============================================================

initParams <- generate_init_lists(4)


# ============================================================
# 7. Run K=4 pleiotropy analysis
# ============================================================

res <- symm_fit_ind_EM_R1(
  t_value = 2,
  testStats = prepped_data$clean_dat,
  initMuList = initParams$initMuList,
  initPiList = initParams$initPiList,
  sameDirAlt = FALSE,
  eps = 10^(-5),
  checkpoint = TRUE
)


# ============================================================
# 8. Process results
# ============================================================

processed <- process_lfdr_results(
  orig_data = prepped_data$extra_dat,
  lfdrResults = res$lfdrResults,
  fdr_threshold = 0.1
)


# ============================================================
# 9. Save results
# ============================================================

write.csv(
  processed$sig_snps_info,
  file.path(output_dir, "csmgmm_pleio_results.csv"),
  row.names = FALSE
)
