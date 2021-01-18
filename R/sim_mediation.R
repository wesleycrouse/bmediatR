#' Simulate Collaborative Cross (CC) data from realized CC genomes
#'
#' This function takes qtl2 formatted genoprobs of CC mice and simulates outcomes/mediators driven by QTL.
#'
#' @param cc_genoprobs A qtl2 formatted genoprobs object for CC mice. 
#' @param cc_map A qtl2 formatted map object that corresponds to the genoprobs.
#' @param num_replicates DEFAULT: 1. The number of strain replicates to use.
#' @param num_sim DEFAULT: 1. The number of simulated outcomes to produce.
#' @param qtl_effect_size DEFAULT = 0.1. The proportion of variation in the simulated outcomes
#' explained by the QTL.
#' @export
#' @examples sim_cc_data()
sim_cc_data <- function(cc_genoprobs, 
                        cc_map,
                        qtl_effect_size = 0.1,
                        strain_effect_size = 0,
                        locus = NULL, 
                        vary_locus = TRUE,
                        num_replicates = 1, 
                        num_sim = 1,
                        M_ID = NULL,
                        sample_method = c("uniform", "crp"),
                        num_alleles = 8, 
                        num_founders = 8,
                        beta = NULL,
                        impute = TRUE,
                        sim_label = "sim_m",
                        return_means = TRUE) {
  
  sample_method <- sample_method[1]
  num_strains <- dim(cc_genoprobs[[1]])[1]
  
  original_effects <- list(qtl_effect_size = qtl_effect_size,
                           strain_effect_size = strain_effect_size)
  
  ############ Scaling effects
  noise_effect_size <- (1 - qtl_effect_size - strain_effect_size)
  reduced_noise_effect_size <- noise_effect_size/num_replicates
  s <- 1/min(c(qtl_effect_size, strain_effect_size, reduced_noise_effect_size)[c(qtl_effect_size, strain_effect_size, reduced_noise_effect_size) != 0])
  qtl_effect_size <- s * qtl_effect_size
  strain_effect_size <- s * strain_effect_size
  noise_effect_size <- s * noise_effect_size
  
  ## Number of individuals
  num_ind <- ifelse(return_means, num_strains, num_strains * num_replicates)
  
  ## Setting number of alleles to match pre-specified QTL effect - convenience
  if (!is.null(beta)) {
    num_alleles <- length(beta)
  }
  
  ## Sampling loci
  loci <- unlist(lapply(cc_genoprobs, function(x) dimnames(x)[[3]]))
  if (is.null(locus)) {
    locus <- sample(loci, size = ifelse(vary_locus, num_sim, 1), replace = TRUE)
    if (vary_locus) {
      locus_index <- 1:num_sim
    }
    else {
      locus_index <- rep(1, num_sim)
    }
  }
  else {
    vary_locus <- FALSE
    locus_index <- rep(1, num_sim)
  }
  
  M <- NULL
  if (!is.null(M_ID)) {
    M <- model_matrix_from_ID(M_ID)
    num_alleles <- length(unique(unlist(strsplit(x = M_ID, split=","))))
  }

  sim_matrix <- matrix(NA, nrow = num_ind, ncol = num_sim)
  for (i in 1:num_sim) {
    this_sim <- sim_cc_qtl(num_replicates = num_replicates,
                           M = M,
                           sample_method = sample_method,
                           num_alleles = num_alleles, 
                           num_founders = num_founders,
                           qtl_effect_size = qtl_effect_size, 
                           strain_effect_size = strain_effect_size,
                           noise_effect_size = noise_effect_size,
                           impute = impute,
                           locus_matrix = qtl2::pull_genoprobpos(genoprobs = cc_genoprobs, marker = locus[locus_index[i]]),
                           return_means = return_means,
                           num_sim = 1,
                           sim_label = sim_label)
    sim_matrix[,i] <- this_sim$data
    if (i == 1) {
      rownames_holder <- rownames(this_sim$data)
    }
  }
  colnames(sim_matrix) <- paste(sim_label, 1:num_sim, sep = "_")
  rownames(sim_matrix) <- rownames_holder
  
  map_df <- qtl2convert::map_list_to_df(cc_map)

  return(list(data = sim_matrix,
              locus = as.character(locus),
              locus_pos = as.numeric(map_df[locus, "pos"]),
              locus_chr = as.character(map_df[locus, "chr"]),
              properties=list(num_alleles = num_alleles,
                              sample_method = sample_method,
                              num_replicates = num_replicates,
                              num_founders = num_founders,
                              qtl_effect_size = original_effects$qtl_effect_size, 
                              strain_effect_size = original_effects$strain_effect_size,
                              impute = impute,
                              M_ID = M_ID,
                              vary_locus = vary_locus,
                              return_means = return_means)))
}

#' Simulate Collaborative Cross (CC) data from a specified locus design matrix
#'
#' This function takes a locus design matrix, i.e., founder haplotypes for CC mice,
#' to simulate outcomes/mediators driven by its QTL.
#'
#' @param locus_matrix A design matrix of founder haplotypes at QTL for CC strains. The number of 
#' rows should correspond to the number of strains. The number of columns is eight for the CC. 
#' @param num_replicates DEFAULT: 1. The number of strain replicates to use.
#' @param num_sim DEFAULT: 1. The number of simulated outcomes to produce.
#' @param qtl_effect_size DEFAULT = 0.1. The proportion of variation in the simulated outcomes
#' explained by the QTL.
#' @export
#' @examples sim_cc_single_locus()
sim_cc_single_locus <- function(locus_matrix, 
                                num_replicates, 
                                num_sim,
                                qtl_effect_size = 0.1,
                                strain_effect_size = 0,
                                M_ID = NULL,
                                sample_method = c("uniform", "crp"),
                                num_alleles = 8, 
                                num_founders = 8,
                                beta = NULL,
                                impute = TRUE,
                                sim_label = "sim_m",
                                return_means = TRUE) {
  
  sample_method <- sample_method[1]
  num_strains <- nrow(locus_matrix)
  
  original_effects <- list(qtl_effect_size = qtl_effect_size,
                           strain_effect_size = strain_effect_size)
  
  ############ Scaling effects
  noise_effect_size <- (1 - qtl_effect_size - strain_effect_size)
  reduced_noise_effect_size <- noise_effect_size/num_replicates
  s <- 1/min(c(qtl_effect_size, strain_effect_size, reduced_noise_effect_size)[c(qtl_effect_size, strain_effect_size, reduced_noise_effect_size) != 0])
  qtl_effect_size <- s * qtl_effect_size
  strain_effect_size <- s * strain_effect_size
  noise_effect_size <- s * noise_effect_size
  
  ## Number of individuals
  num_ind <- ifelse(return_means, num_strains, num_strains * num_replicates)
  
  ## Setting number of alleles to match pre-specified QTL effect - convenience
  if (!is.null(beta)) {
    num_alleles <- length(beta)
  }
  
  M <- NULL
  if (!is.null(M_ID)) {
    M <- model_matrix_from_ID(M_ID)
    num_alleles <- length(unique(unlist(strsplit(x = M_ID, split=","))))
  }
  
  sim_matrix <- matrix(NA, nrow = num_ind, ncol = num_sim)
  for (i in 1:num_sim) {
    this_sim <- sim_cc_qtl(num_replicates = num_replicates,
                           M = M,
                           beta = beta,
                           sample_method = sample_method,
                           num_alleles = num_alleles, 
                           num_founders = num_founders,
                           qtl_effect_size = qtl_effect_size, 
                           strain_effect_size = strain_effect_size,
                           noise_effect_size = noise_effect_size,
                           impute = impute,
                           locus_matrix = locus_matrix,
                           return_means = return_means,
                           num_sim = 1,
                           sim_label = sim_label)
    sim_matrix[,i] <- this_sim$data
    if (i == 1) {
      rownames_holder <- rownames(this_sim$data)
    }
  }
  colnames(sim_matrix) <- paste(sim_label, 1:num_sim, sep = "_")
  rownames(sim_matrix) <- rownames_holder
  
  return(list(data = sim_matrix,
              locus_matrix = locus_matrix,
              properties=list(num_alleles = num_alleles,
                              sample_method = sample_method,
                              num_replicates = num_replicates,
                              num_founders = num_founders,
                              qtl_effect_size = original_effects$qtl_effect_size, 
                              strain_effect_size = original_effects$strain_effect_size,
                              impute = impute,
                              M_ID = M_ID,
                              return_means = return_means)))
}

sim_cc_qtl <- function(locus_matrix,
                       num_replicates,
                       M = NULL,
                       sample_method = c("uniform", "crp"),
                       qtl_effect_size, 
                       beta = NULL,
                       strain_effect_size,
                       noise_effect_size,
                       num_alleles = 8, 
                       num_founders = 8,
                       num_sim,
                       impute = TRUE, 
                       return_means = TRUE,
                       sim_label = "sim_m",
                       ...) {
  
  sample_method <- sample_method[1]
  strains <- rownames(locus_matrix)
  
  ## Imputation for simulation
  if (impute) {
    locus_matrix <- t(apply(locus_matrix, 1, function(x) rmultinom(1, 1, x)))
    rownames(locus_matrix) <- strains
  }
  D <- locus_matrix ## Saved for variance calculations
  ZD <- D[rep(1:nrow(locus_matrix), each = num_replicates),]
  
  ## Strain
  Z <- incidence_matrix(factor(strains, levels = strains))
  Z <- Z[rep(1:nrow(Z), each = num_replicates),]
  if (strain_effect_size != 0) {
    strain_effect <- rnorm(n = nrow(locus_matrix))
    strain_effect <- (strain_effect - mean(strain_effect))/sqrt(non_sample_var(strain_effect))
    strain_effect <- strain_effect * sqrt(strain_effect_size)
    
    strain_predictor <- Z %*% matrix(strain_effect, ncol = 1)
  }
  else { strain_predictor <- rep(0, nrow(Z)) }
  
  ## QTL
  if (qtl_effect_size != 0) {
    qtl_effect <- sim_qtl_model_and_effects(num_alleles = num_alleles, 
                                            num_founders = num_founders, 
                                            M = M,
                                            sample_method = sample_method,
                                            beta = beta,
                                            ...)
  }
  else {
    qtl_effect <- list(M = diag(8),
                       beta = rep(0, 8))
  }
  M <- qtl_effect$M
  raw_beta <- qtl_effect$beta 
  if (qtl_effect_size != 0) {
    beta <- (raw_beta - mean(raw_beta))/sqrt(non_sample_var(raw_beta))
  }
  else { 
    beta <- raw_beta 
  }
  
  ## Scaling (ZMB)
  var_ratio <- ifelse(qtl_effect_size != 0,
                      c(non_sample_var(Z %*% D %*% M %*% beta)/non_sample_var(beta)),
                      1)
  
  if (var_ratio != 0) { # Case when more than one allele is observed
    beta <- beta*sqrt(qtl_effect_size)*sqrt(1/var_ratio)
  }
  else { # Case when only one allele is observed, should result in a null scan
    beta <- beta*sqrt(qtl_effect_size)
  }
  
  scaled_qtl_effects <- calc_qtl_effect(beta = beta,
                                        M = M,
                                        D = D,
                                        Z = Z,
                                        strain_effect_size = strain_effect_size,
                                        noise_effect_size = noise_effect_size)

  qtl_predictor <- scaled_qtl_effects$qtl_predictor
  
  sim_data <- matrix(NA, nrow = nrow(Z), ncol = num_sim)
  for (i in 1:num_sim) {
    scaled_resid <- calc_scaled_residual(noise_effect_size = noise_effect_size,
                                         n = nrow(Z))
    sim_data[,i] <- qtl_predictor + strain_predictor + scaled_resid
  }

  colnames(sim_data) <- paste(sim_label, 1:ncol(sim_data), sep = "_")
  rownames(sim_data) <- paste(rep(strains, each = num_replicates), 1:num_replicates, sep = "_")
  
  if (return_means) {
    sim_data <- data.frame(strain = rep(strains, each = num_replicates), sim_data) %>%
      tidyr::gather(key = "sim", value = "phenotype", -strain) %>%
      dplyr::group_by(strain, sim) %>%
      dplyr::summarize(phenotype = mean(phenotype)) %>%
      dplyr::ungroup() %>%
      tidyr::spread(key = "sim", value = "phenotype") %>%
      tibble::column_to_rownames("strain") %>%
      as.matrix
  }
  return(list(data = sim_data,
              properties = list(qtl_effect_size = qtl_effect_size,
                                strain_effect_size = strain_effect_size,
                                num_alleles = num_alleles,
                                sample_method = sample_method,
                                num_replicates = num_replicates,
                                impute = impute,
                                return_means = return_means)))
}

## From Wes, returns SDP matrix and QTL effects
sim_qtl_model_and_effects <- function(num_alleles = 8, 
                                      num_founders = 8, 
                                      M = NULL,
                                      sample_method = c("uniform", "crp"), 
                                      beta = NULL,
                                      ...){
  
  sample_method <- sample_method[1]
  if (is.null(M)) {
    M <- matrix(0, num_founders, num_alleles)
    M[cbind(sample(1:num_founders, num_alleles), 1:num_alleles)] <- 1
    
    if (sample_method == "uniform") {
      M[which(rowSums(M)==0),] <- t(rmultinom(num_founders - num_alleles, 1, rep(1/num_alleles, num_alleles)))
    }
    else if (sample_method == "crp") {
      for (i in which(rowSums(M)==0)) {
        M[i,] <- rmultinom(1, 1, colSums(M)/sum(M))
      }
    }
  }
  
  if (is.null(beta)) {
    beta <- rnorm(num_alleles)
  }
  
  effect <- list(M = M, 
                 beta = beta)
  effect
}


#' Produce SDP mapping matrix from SDP string
#'
#' This function takes an SDP string, i.e., "0,0,0,0,1,1,1,1", and converts it into the a
#' matrix that maps the eight allele haplotype design matrix to < eight functional allele
#' design matrix.
#'
#' @param M_ID SDP string. Example is "0,0,0,0,1,1,1,1", which represent two functional 
#' alleles evenly distributed across the eight founder strains.
#' 
#' @export
#' @examples model_matrix_from_ID()
model_matrix_from_ID <- function(M_ID) {
  
  m <- as.numeric(unlist(strsplit(M_ID, ","))) + 1
  J <- length(m)
  M <- matrix(0, J, max(m))
  M[cbind(1:J, m)] <- 1
  M
}

## Produce an incidence matrix from a factor
incidence_matrix <- function(fact) {
  
  m <- diag(nlevels(fact))[fact,]
  colnames(m) <- levels(fact)
  return(m)
}

## From miqtl originally
straineff_mapping_matrix <- function(M = 8) {
  
  T <- M*(M+1)/2
  mapping <- matrix(rep(0, T*M), M, T)
  idx <- 1;
  for (i in 1:M){
    mapping[i, idx] <- mapping[i, idx] + 2
    idx <- idx + 1;
  }
  for (i in 2:M){
    for (j in 1:(i-1)){
      mapping[i, idx] <- mapping[i, idx] + 1;
      mapping[j, idx] <- mapping[j, idx] + 1;
      idx <- idx + 1;
    }
  }
  t(mapping)
}

non_sample_var <- function(x) {
  
  var_x <- var(x) * ((length(x) - 1)/length(x))
  var_x
}

calc_qtl_effect <- function(beta,
                            M,
                            D,
                            Z,
                            strain_effect_size,
                            noise_effect_size) {
  
  MB <- M %*% beta
  DMB <- D %*% MB
  ZDMB <- Z %*% DMB
  
  B_effect <- non_sample_var(beta)
  MB_effect <- non_sample_var(MB)
  DMB_effect <- non_sample_var(DMB)
  ZDMB_effect <- non_sample_var(ZDMB)
  
  summaries <- c(B_effect, MB_effect, DMB_effect, ZDMB_effect,
                 B_effect/(B_effect + strain_effect_size + noise_effect_size),
                 MB_effect/(MB_effect + strain_effect_size + noise_effect_size),
                 DMB_effect/(DMB_effect + strain_effect_size + noise_effect_size),
                 ZDMB_effect/(ZDMB_effect + strain_effect_size + noise_effect_size))
  names(summaries) <- c("B_effect", "MB_effect", "DMB_effect", "ZDMB_effect",
                        "B_ve", "MB_ve", "DMB_ve", "ZDMB_ve")
  
  return(list(qtl_predictor = ZDMB,
              summaries = summaries))
}

## Draws and scales residuals in single function
calc_scaled_residual <- function(noise_effect_size, n) {
  
  residual <- rnorm(n = n)
  residual <- (residual - mean(residual))/sqrt(non_sample_var(residual))
  residual <- residual*sqrt(noise_effect_size)
  residual
}

#' Simulate target from mediator
#'
#' This function takes a mediator matrix and simulates downstream
#' outcomes or targets.
#'
#' @param M Vector or matrix of mediator variables. Multiple mediator variables are supported. 
#' Names or rownames must match across y and X, and Z and w (if provided).
#' @param mediator_effect_size DEFAULT = 0.5. The proportion of variation in the simulated outcomes
#' explained by the mediator.
#' 
#' @export
#' @examples sim_target()
sim_target <- function(M,
                       mediator_effect_size = 0.5,
                       target_to_mediator_var_ratio = 1,
                       sim_label = "sim_y") {
  
  Y <- matrix(0, nrow = nrow(M), ncol = ncol(M))
  rownames(Y) <- rownames(M)
  colnames(Y) <- paste(sim_label, 1:ncol(M), sep = "_")
  for (i in 1:ncol(M)) {
    m_var <- non_sample_var(M[,i])
    y_var <- m_var * (1 - mediator_effect_size)
    e <- rnorm(n = nrow(M), mean = 0, sd = 1)
    if (mediator_effect_size != 0) {
      e_var <- ((1 - mediator_effect_size) * m_var)/mediator_effect_size
      Y[,i] <- M[,i] * sqrt(target_to_mediator_var_ratio)
    } else{
      e_var <- m_var
    }
    
    e <- (e * sqrt(e_var))/sqrt(non_sample_var(e))
    Y[,i] <- Y[,i] + e * sqrt(target_to_mediator_var_ratio)
  }
  Y
}

#' Simulate a balanced eight allele design matrix
#'
#' This function simulates a balanced eight allele design 
#' matrix for a specified number of observations of each
#' founder allele.
#'
#' @param founder_allele_reps DEFAULT: 10. The number of observations of each founder allele.
#' For example, the default of 10 will produce a design matrix representing 80 artificial CC
#' strains.
#' 
#' @export
#' @examples sim_balanced_locus()
sim_balanced_locus <- function(founder_allele_reps = 10) {
  
  num_ind <- 8 * founder_allele_reps
  strain_num <- sapply(1:num_ind, function(i) ifelse(nchar(i) == 1, paste0("00", i), i))
  strain_num <- sapply(1:num_ind, function(i) ifelse(nchar(strain_num[i]) == 2, paste0("0", i), strain_num[i]))
  
  locus_matrix <- matrix(0, nrow = num_ind, ncol = 8)
  colnames(locus_matrix) <- LETTERS[1:8]
  rownames(locus_matrix) <- paste0("CC", strain_num)
  
  hap_sets <- split(sample(1:num_ind), ceiling(seq_along(1:num_ind)/founder_allele_reps))
  
  for (i in 1:length(hap_sets)) {
    locus_matrix[hap_sets[[i]],i] <- 1
  }
  
  locus_matrix
}



