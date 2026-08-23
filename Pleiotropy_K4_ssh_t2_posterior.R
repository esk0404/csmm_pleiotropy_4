library(dplyr)
library(mvtnorm)
library(data.table)
library(bindata)
library(magrittr)
library(devtools)
library(rje)
library(ks)
library(csmGmm)
library(here)
library(locfdr)
library(stats4)
library(BiocGenerics)
library(generics)
library(purrr)
library(readr)
library(GenomicRanges)
library(rtracklayer)
library(R.utils)

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


calc_dens_ind_multiple <- function(x, Zmat) {
  K <- ncol(Zmat)
  tempSum <- rep(0, nrow(Zmat))
  for (k_it in 1:K) {
    tempSum <- tempSum + stats::dnorm(Zmat[, k_it], mean=x[k_it], sd=1, log=TRUE)
  }
  exp(tempSum)
}



find_max_means_R1 <- function(muInfo) {
  
  # iterate, skip the first (0) and last (alternative)
  listLength <- length(muInfo)
  K <- nrow(muInfo[[1]])
  #S1 <- c(1,2,4,8)+1  # t=1
  S1 <- c(1:6,8:10,12)+1   # t=2
  # Com_Set <- compute_sets(K,t)
  # S1 <- Com_Set[1]
  # just keep finding the max
  maxMeans <- rep(0, K) 
  for (element_it in S1) {
    tempMat <- cbind(muInfo[[element_it]], maxMeans)
    maxMeans <- apply(tempMat, 1, max)
  }
  # return K*1 vector
  return(maxMeans)
}



symm_fit_ind_EM_R1 <- function(t_value, testStats, initMuList, initPiList, sameDirAlt=FALSE, eps = 10^(-5), checkpoint=TRUE) {
  
  # number of composite null hypotheses
  J <- nrow(testStats)
  # number of dimensions
  K <- ncol(testStats)
  B <- 2^K - 1
  # number of hl configurations - 1
  L <- 3^K - 1
  #number of the non_zero values
  t <- t_value
  # make all configurations
  Hmat <- expand.grid(rep(list(-1:1), K))
  # attach the bl
  blVec <- rep(0, nrow(Hmat))
  slVec <- rep(0, nrow(Hmat))
  for (k_it in K:1) {
    blVec <- blVec + 2^(K - k_it) * abs(Hmat[, k_it])
    slVec <- slVec + abs(Hmat[, k_it]) #number of non_zeros
  }
  # symmetric alternative
  sum_Hmat <- apply(Hmat[, 1:K], 1, sum)
  symAltVec <- ifelse((t+1 <= sum_Hmat & sum_Hmat <= K) |  (sum_Hmat >= -K & sum_Hmat <= -t-1), 1, 0)
  #symAltVec <- ifelse(apply(Hmat[, 1:K], 1, sum) == K | apply(Hmat[, 1:K], 1, sum) == -K, 1, 0)
  
  # sort Hmat
  Hmat <- Hmat %>% dplyr::mutate(bl = blVec) %>%
    dplyr::mutate(sl = slVec) %>%
    dplyr::mutate(symAlt = symAltVec) %>%
    dplyr::arrange(.data$bl, across(starts_with("Var"))) %>%
    dplyr::mutate(l = 0:(nrow(.) - 1)) %>%
    dplyr::relocate(.data$l, .before = .data$symAlt)
  
  # initialize
  # the muInfo and piInfo are lists that hold the information in a more compact manner.
  # the allMu and allPi are matrices that repeat the information so the calculations can be performed faster.
  muInfo <- initMuList
  piInfo <- initPiList
  oldParams <- c(unlist(piInfo), unlist(muInfo))
  MbVec <- sapply(piInfo, FUN=length)
  
  # run until convergence
  diffParams <- 10
  iter <- 0
  while (diffParams > eps) {
    
    #####################################################
    # first, update allPi and allMu with the new parameter values
    
    # allPi holds the probabilities of each configuration (c 1), the bl of that
    # configuration (c 2), the sl of that configuration (c3), the m of that configuration (c 4),
    # the l of that configuration (c 5), and whether it's part of the symmetric alternative (c 6).
    
    # number of rows in piMat is L if Mb = 1 for all b.
    # number of rows is \sum(l = 0 to L-1) {Mbl}
    allPi <- c(piInfo[[1]], 0, 0, 1, 0, 0)
    
    # allMu holds the mean vectors for each configuration in each row, number of columns is K
    allMu <- rep(0, K)
    
    # loop through possible values of bl
    for (b_it in 1:B) {
      
      # Hmat rows with this bl
      tempH <- Hmat %>% dplyr::filter(.data$bl == b_it)
      
      # loop through possible m for this value of bl
      for (m_it in 1:MbVec[b_it + 1]) {
        allPi <- rbind(allPi, cbind(rep(piInfo[[b_it + 1]][m_it] / (2^tempH$sl[1]), nrow(tempH)), tempH$bl,
                                    tempH$sl, rep(m_it, nrow(tempH)), tempH$l, tempH$symAlt))
        for (h_it in 1:nrow(tempH)) {
          allMu <- cbind(allMu, unlist(tempH %>% dplyr::select(-.data$bl, -.data$sl, -.data$l, -.data$symAlt) %>%
                                         dplyr::slice(h_it)) * muInfo[[b_it + 1]][, m_it])
        }
      } # done looping through different m
    } # dont looping through bl
    colnames(allPi) <- c("Prob", "bl", "sl", "m", "l", "symAlt")
    
    ##############################################################################################
    # this is the E step where we calculate Pr(Z|c_l,m) for all c_l,m.
    # each l,m is one column.
    # only independence case for now
    if (ncol(testStats) == 2) {
      conditionalMat <- sapply(X=data.frame(allMu), FUN=calc_dens_ind_2d, Zmat = testStats) %>%
        sweep(., MARGIN=2, STATS=allPi[, 1], FUN="*")
    } else if (ncol(testStats) == 3) {
      conditionalMat <- sapply(X=data.frame(allMu), FUN=calc_dens_ind_3d, Zmat = testStats) %>%
        sweep(., MARGIN=2, STATS=allPi[, 1], FUN="*")
    } else {
      conditionalMat <- sapply(X=data.frame(allMu), FUN=calc_dens_ind_multiple, Zmat = testStats) %>%
        sweep(., MARGIN=2, STATS=allPi[, 1], FUN="*")
    }
    probZ <- apply(conditionalMat, 1, sum)
    AikMat <- conditionalMat %>% sweep(., MARGIN=1, STATS=probZ, FUN="/")
    Aik_alln <- apply(AikMat, 2, sum) / J
    
    ###############################################################################################
    # this is the M step for probabilities of hypothesis space
    for (b_it in 0:B) {
      for (m_it in 1:MbVec[b_it + 1]) {
        tempIdx <- which(allPi[, 2] == b_it & allPi[, 4] == m_it)
        piInfo[[b_it + 1]][m_it] <- sum(Aik_alln[tempIdx])
      }
    }
    
    # M step for the means
    # loop through values of bl
    # do the alternative last, enforce that it must be larger in magnitude than the nulls
    for (b_it in 1:B) {
      
      tempHmat <- Hmat %>% dplyr::filter(.data$bl == b_it)
      # loop through m
      for (m_it in 1:MbVec[b_it + 1]) {
        tempMuSum <- rep(0, nrow(allMu))
        tempDenom <- rep(0, nrow(allMu))
        
        # these are the classes that contribute to \mu_bl,m
        AikIdx <- which(allPi[, 2] == b_it & allPi[, 4] == m_it)
        for (idx_it in 1:length(AikIdx)) {
          tempAik <- AikIdx[idx_it]
          tempHvec <- tempHmat %>% dplyr::select(-.data$bl, -.data$sl, -.data$l, -.data$symAlt) %>%
            dplyr::slice(idx_it) %>% unlist(.)
          
          tempMuSum <- tempMuSum + colSums(AikMat[, tempAik] * sweep(x = testStats, MARGIN = 2,
                                                                     STATS = tempHvec, FUN="*"))
          tempDenom <- tempDenom + rep(J * Aik_alln[tempAik], length(tempDenom)) * abs(tempHvec)
        } # done looping for one l, m
        whichZero <- which(tempDenom == 0)
        tempDenom[whichZero] <- 1
        muInfo[[b_it + 1]][, m_it] <- tempMuSum / tempDenom
        
        # make sure mean constraint is satisfied
        # S2 <- c(3,5,6,7,9,10,11,12,13,14,15)  # t=1
        S2 <- c(7,11,13,14,15)  # t=2
        #Com_Set <- compute_sets(K,t)
        #S2 <- Com_Set[2]
        
        if (b_it %in% S2) {
          maxMeans <- find_max_means_R1(muInfo)
          whichSmaller <- which(muInfo[[b_it+1]][, m_it] < maxMeans)
          if (length(whichSmaller) > 0) {
            muInfo[[b_it+1]][whichSmaller, m_it] <- maxMeans[whichSmaller]
            #if (b_it == B) {
            
            #maxMeans <- find_max_means_R1(muInfo)
            #whichSmaller <- which(muInfo[[b_it + 1]][, m_it] < maxMeans)
            #if (length(whichSmaller) > 0) {
            #  muInfo[[b_it + 1]][whichSmaller, m_it] <- maxMeans[whichSmaller]
            
          }
        } # done with mean constraint
        
      } # done looping through m
      
    } # done updating means
    
    ###############################################################################################
    # find difference
    allParams <- c(unlist(piInfo), unlist(muInfo))
    diffParams <- sum((allParams - oldParams)^2)
    
    # update
    oldParams <- allParams
    iter <- iter + 1
    if (checkpoint) {
      cat(iter, " - ", diffParams, "\n", allParams, "\n")
    }
  }
  
  # calculate local fdrs
  if (sameDirAlt) {
    nullCols <- which(allPi[, 6] == 0)
  } else {
    nullCols <- which(allPi[, 3] <= t)  # revised on 2025/09/28
    #nullCols <- which(allPi[, 3] < K)
  }
  probNull <- apply(conditionalMat[, nullCols], 1, sum)
  lfdrResults <- probNull / probZ
  
  return(list(piInfo = piInfo, muInfo = muInfo, iter = iter,
              lfdrResults = lfdrResults,
              conditionalMat = conditionalMat,
              AikMat = AikMat,
              allPi = allPi,
              allMu = allMu))
}

#------------------------------------------------------------------------------------------
here::i_am("Pleiotropy_K4_ssh_t2_posterior.R")

# define paths
dataDir <- here::here("Data")
outputDir <- here::here("output")

# liftover chain
chain_38to37 <- import.chain(
  file.path(dataDir, "hg38ToHg19.over.chain")
)


# read data
prostate_read <- fread(file.path(dataDir, "GCST006085.txt"))
prostate <- prostate_read %>% select(Chr, position, rs_id, SNP,
                                     Allele2, Allele1,
                                     Effect, StdErr, Pvalue)

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

breast_read <- fread(file.path(dataDir, "GCST004988.txt"))
breast <- breast_read %>% select(var_name,chr, 
                                 phase3_1kg_id,
                                 position_b37,
                                 bcac_gwas_all_beta, 
                                 bcac_gwas_all_P1df,
                                 bcac_gwas_all_se, a0, a1)

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

lung_read <- fread(file.path(dataDir, "GCST004744_buildGRCh37.tsv"))
lung <- lung_read %>% select(variant_id, chromosome,
                             base_pair_location,
                             beta,
                             standard_error,
                             effect_allele,
                             other_allele,
                             p_value)

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


ovarian_read <- fread(file.path(dataDir, "GCST90455658.tsv.gz"))
ovarian <- ovarian_read %>% select(variant_id, chromosome, base_pair_location,
                                   other_allele, effect_allele,
                                   beta, standard_error, p_value) # build 38



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


rm(chain_38to37, lifted, df, gr, lifted_gr, new_chr, keep0, keep1)

datasets <- list(breast_new,lung_new,ovarian_liftover,prostate_new)
prepped_data <- prepare_csmgmm_data(datasets=datasets, dataset_names = c("Breast", "Lung", "Ovarian", "Prostate")) 


nDims <- 4
nComp <- 2^nDims
initPiList <- vector("list", nComp)
initPiList[[1]] <- 0.90

for (i in 2:nComp) {initPiList[[i]] <- 0.10 / (nComp - 1)}

# the symm_fit_ind.R code will add the appropriate 0s to initMuList
initMuList <- list(matrix(data=rep(0, nDims), nrow=nDims, ncol=1))
for (i in 2:(2^nDims)) { 
  initMuList[[i]] <- matrix(data=rep(4, nDims), nrow=nDims, ncol=1) 
}                                               


res <- symm_fit_ind_EM_R1(t_value = 2, 
                          testStats = prepped_data$clean_dat, 
                          initMuList = initMuList, 
                          initPiList = initPiList, 
                          sameDirAlt = FALSE,
                          eps = 10^(-5), checkpoint = TRUE)



conditionalMat <- res$conditionalMat
AikMat <- res$AikMat
prior <- res$allPi[, "Prob"]


conditional_df <- as.data.frame(conditionalMat)
colnames(conditional_df) <- paste0(
  "conditional_", seq_len(ncol(conditionalMat))
)

posterior_df <- as.data.frame(AikMat)
colnames(posterior_df) <- paste0(
  "posterior_", seq_len(ncol(AikMat))
)


prior_df <- as.data.frame(t(prior))
colnames(prior_df) <- paste0("prior_", seq_along(prior))

# --------------------------------------------------
# Calculate cumulative average lfdr
# --------------------------------------------------

sorted_idx <- order(res$lfdrResults)

cum_avg <- cumsum(
  res$lfdrResults[sorted_idx]
) / seq_along(sorted_idx)

cum_avg_original <- numeric(length(res$lfdrResults))

cum_avg_original[sorted_idx] <- cum_avg

# --------------------------------------------------
# Add lfdr and cumulative average
# --------------------------------------------------

all_results <- prepped_data$extra_dat %>%
  mutate(
    lfdrResults = res$lfdrResults,
    cum_avg = cum_avg_original
  ) %>%
  bind_cols(
    conditional_df,
    posterior_df
  )

# --------------------------------------------------
# 7. Select significant SNPs
# --------------------------------------------------

sig_snps_info <- all_results %>%
  filter(cum_avg < 0.1) %>%
  arrange(lfdrResults) %>%
  select(Chr, BP, EA, OA, 
    rsid,
    lfdrResults,
    cum_avg,
    starts_with("conditional_"),
    starts_with("posterior_")
  )

# --------------------------------------------------
# 8. Write output
# --------------------------------------------------

posterior_outfile <- file.path(
  outputDir,
  "Pleiotropy_K4_ssh_posteriors.csv"
)

prior_outfile <- file.path(
  outputDir,
  "Pleiotropy_K4_ssh_priors.csv"
)

write.csv(
  sig_snps_info,
  posterior_outfile,
  row.names = FALSE
)

write.csv(
  prior_df,
  prior_outfile,
  row.names = FALSE
)