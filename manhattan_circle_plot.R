library(readr)
library(dplyr)
library(rtracklayer)
library(GenomicRanges)
library(data.table)
library(circlize)
library(scales)

# ----------------------------
# 0. Load data
# ----------------------------
setwd("/Users/ekim14/Downloads")

prostate <- read_tsv("GCST006085.txt",
                     col_select = c(Chr, position, Pvalue),
                     progress = FALSE)

breast <- read_tsv("GCST004988.txt",
                   col_select = c(chr, 
                                  position_b37,
                                  bcac_gwas_all_P1df),
                   progress = FALSE)

lung <- read_tsv("GCST004744_buildGRCh37.tsv",
                 col_select = c(chromosome,
                                base_pair_location,
                                p_value))

ovarian <- read_tsv("GCST90455658.tsv.gz",
                    col_select = c(chromosome, base_pair_location, p_value))

# ----------------------------
# 1. Preprocess datasets
# ----------------------------
lung_new <- lung %>%
  transmute(chr = as.factor(chromosome),
            logp = -log10(p_value),
            position = base_pair_location) %>%
  filter(!is.na(logp), !is.na(chr), !is.na(position), chr != 23) %>%
  mutate(chr = droplevels(chr))
            
rm(lung)

breast_new <- breast %>% 
  transmute(chr = as.factor(chr),
            p_value = as.numeric(bcac_gwas_all_P1df),
            position = position_b37) %>%
  mutate(logp = -log10(p_value)) %>%
  filter(!is.na(logp), !is.na(chr), !is.na(position), chr != 23) %>%
  mutate(chr = droplevels(chr))

rm(breast)

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


rm(chain_38to37)

ovarian_new <- ovarian_liftover %>%
  transmute(chr = as.factor(chromosome),
            position = as.numeric(base_pair_location),
            logp = -log10(p_value)) %>%
  filter(!is.na(logp), !is.na(chr), !is.na(position), chr != 23) %>%
  mutate(chr = droplevels(chr))

rm(lifted, df, gr, ovarian_liftover, lifted_gr, new_chr, keep0, keep1)


prostate_new <- prostate %>%
  transmute(chr = as.factor(Chr),
         position = position,
         logp = -log10(Pvalue)) %>%
  filter(!is.na(logp), !is.na(chr), !is.na(position), chr != 23) %>%
  mutate(chr = droplevels(chr))

rm(prostate)

# csmGmm dataset
getwd()
setwd("/Users/ekim14/Downloads")
pleio <- fread("~/Downloads/Pleiotropy_all_chr_plot.csv")


pleio_new <- pleio %>%
  transmute(
    chr = as.factor(Chr),
    position = Pos,
    loglfdr = -log10(lfdrResults)
  ) %>%
  filter(!is.na(chr), !is.na(position), !is.na(loglfdr))


# list datasets

datasets <- list(prostate_new, ovarian_new, breast_new, lung_new)


bin_dataset <- function(df, bin_size = 50000) {
  df %>%
    mutate(bin = position %/% bin_size) %>%
    group_by(chr, bin) %>%
    slice_max(logp, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(position = bin * bin_size)
}

datasets_binned <- lapply(datasets, bin_dataset)



# get real chromosome lengths
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
library(BSgenome.Hsapiens.UCSC.hg19)

chr_lengths <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)
chromosomes <- paste0("chr", c(1:22))
#chromosomes <- paste0("chr", c(1:22, "X"))
chr_lengths <- chr_lengths[chromosomes]
#names(chr_lengths)[names(chr_lengths) == "chrX"] <- "chr23"
names(chr_lengths) <- gsub("^chr", "", names(chr_lengths))
chr_lengths

# Start plotting
circos.clear()
circos.par(
  start.degree = 90,
  track.margin = c(0.006, 0.006)
)


circos.initialize(
  factors = names(chr_lengths),
  xlim = cbind(rep(0, length(chr_lengths)), chr_lengths)
)

circos.trackPlotRegion(
  ylim = c(0, 1),
  track.height = 0.09,
  bg.border = NA,
  panel.fun = function(x, y) {
    
    chr <- CELL_META$sector.index
    
    circos.text(
      CELL_META$xcenter,
      0.08,   #  push outward
      chr,
      facing = "bending.outside",
      niceFacing = TRUE,
      cex = 0.7
    )
  }
)



track_colors <- c("#A7B8D6", "#9DC3A6", "#E6B89C", "#C8A2C8")
names(datasets_binned) <- c("Prostate", "Ovarian", "Breast", "Lung")

global_ymax <- max(
  sapply(datasets_binned, function(df) max(df$logp, na.rm = TRUE))
)

global_ymax <- ceiling(global_ymax)
ticks <- c(0, global_ymax/2, global_ymax)


for (i in seq_along(datasets_binned)) {
  
  df <- datasets_binned[[i]]
  col_i <- track_colors[i]
  
  circos.trackPlotRegion(
    factors = df$chr,
    y = df$logp,
    ylim = c(0, global_ymax),
    track.height = 0.12,
    bg.border = NA,
    
    panel.fun = function(x, y) {
      
      chr <- CELL_META$sector.index
      idx <- df$chr == chr
      
      circos.points(
        df$position[idx],
        df$logp[idx],
        col = col_i,
        pch = 16,
        cex = 0.5
      )
      
      if (CELL_META$sector.numeric.index == 1) {
        circos.yaxis(
          side = "left",
          at = ticks,
          labels = ticks,
          labels.cex = 0.5,
          tick.length = convert_y(2, "mm"),
          labels.niceFacing = TRUE
        )
      }
    }
  )
}

ymax_lfdr <- ceiling(max(pleio_new$loglfdr, na.rm = TRUE))
ticks_lfdr <- c(0, ymax_lfdr/2, ymax_lfdr)

circos.trackPlotRegion(
  ylim = c(0, ymax_lfdr),
  track.height = 0.2,
  bg.border = NA,
  
  panel.fun = function(x, y) {
    
    chr <- CELL_META$sector.index
    idx <- pleio_new$chr == chr
    
    circos.points(
      pleio_new$position[idx],
      pleio_new$loglfdr[idx],
      col = "red",
      pch = 16,
      cex = 0.7
    )
    
    if (CELL_META$sector.numeric.index == 1) {
      circos.yaxis(
        side = "left",
        at = ticks_lfdr,
        labels = ticks_lfdr,
        labels.cex = 0.5,
        tick.length = convert_y(3, "mm"),
        labels.niceFacing = TRUE
      )
    }
  }
)

circos.trackPlotRegion(
  ylim = c(0, 1),
  track.height = 0.08,
  bg.border = NA,
  panel.fun = function(x, y) {
    
    circos.rect(
      CELL_META$xlim[1], 0,
      CELL_META$xlim[2], 1,
      col = "pink",
      border = NA
    )
  }
)

legend(
  "topright",
  legend = c(names(datasets_binned), "Pleiotropic Effects"),
  col = c(track_colors, "red"),
  pch = 16,
  pt.cex = 1,
  cex = 0.7,
  bty = "l"
)

