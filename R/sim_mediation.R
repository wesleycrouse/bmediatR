#' Simulate target (y) and mediator (m) for partial mediation from a specified locus design matrix
#'
#' This function takes a locus design matrix, e.g., founder haplotypes dosages from CC mice, 
#' and simulates mediators and downstream outcomes or targets.
#'
#' @param locus_matrix A design matrix of founder haplotypes at QTL for MPP founder strains. The number of 
#' rows should correspond to the number of strains. The number of columns is eight for the CC/DO. 
#' @param num_replicates DEFAULT: 1. The number of strain replicates to use.
#' @param num_sim DEFAULT: 1. The number of simulated outcomes to produce.
#' @param qtl_a_effect_size DEFAULT = 0.1. The proportion of variation in the simulated mediators
#' explained by the QTL.
#' @param m_b_effect_size DEFAULT = 0.5. The proportion of variation in the simulated outcomes
#' explained by the mediator.
#' @param qtl_c_effect_size DEFAULT = 0.1. The proportion of variation in the simulated outcomes
#' explained by the direct effect of the QTL.
#' 
#' @export
#' @examples sim_mpp_single_locus_partial()
sim_mpp_single_locus_partial <- function(locus_matrix, 
                                         num_replicates, 
                                         num_sim,
                                         qtl_a_effect_size = 0.1,
                                         m_b_effect_size = 0.5,
                                         qtl_c_effect_size = 0.1,
                                         beta_a = NULL,
                                         beta_c = NULL,
                                         M_ID_a = NULL,
                                         M_ID_c = NULL,
                                         num_alleles_a = 8, 
                                         num_alleles_c = 8,
                                         m_strain_effect_size = 0,
                                         y_strain_effect_size = 0,
                                         sample_method = c("uniform", "crp"),
                                         num_founders = 8,
                                         impute = TRUE,
                                         sim_m_label = "sim_m",
                                         sim_y_label = "sim_y",
                                         return_means = TRUE) {
  
  sample_method <- sample_method[1]
  num_strains <- nrow(locus_matrix)
  
  original_effects <- list(qtl_a_effect_size = qtl_a_effect_size,
                           m_b_effect_size = m_b_effect_size,
                           qtl_c_effect_size = qtl_c_effect_size,
                           m_strain_effect_size = m_strain_effect_size,
                           y_strain_effect_size = y_strain_effect_size)
  
  ############ Scaling effects
  m_noise_effect_size <- (1 - qtl_a_effect_size - m_strain_effect_size)
  m_reduced_noise_effect_size <- m_noise_effect_size/num_replicates
  m_scale_factor <- 1/sum(c(qtl_a_effect_size, m_strain_effect_size, m_reduced_noise_effect_size))
  m_qtl_a_effect_size <- m_scale_factor * qtl_a_effect_size
  m_strain_effect_size <- m_scale_factor * m_strain_effect_size
  m_noise_effect_size <- m_scale_factor * m_noise_effect_size
  
  y_noise_effect_size <- (1 - m_b_effect_size - qtl_c_effect_size - y_strain_effect_size)
  y_reduced_noise_effect_size <- y_noise_effect_size/num_replicates
  y_scale_factor <- 1/sum(c(m_b_effect_size, qtl_c_effect_size, y_strain_effect_size, y_reduced_noise_effect_size))
  y_m_b_effect_size <- y_scale_factor * m_b_effect_size
  y_qtl_c_effect_size <- y_scale_factor * qtl_c_effect_size
  y_strain_effect_size <- y_scale_factor * y_strain_effect_size
  y_noise_effect_size <- y_scale_factor * y_noise_effect_size
  
  ## Number of individuals
  num_ind <- ifelse(return_means, num_strains, num_strains * num_replicates)
  
  ## Setting number of alleles to match pre-specified QTL effect - convenience
  if (!is.null(beta_a)) {
    num_alleles_a <- length(beta_a)
  }
  if (!is.null(beta_c)) {
    num_alleles_c <- length(beta_c)
  }
  
  M_a <- NULL
  if (!is.null(M_ID_a)) {
    M_a <- model_matrix_from_ID(M_ID_a)
    num_alleles_a <- length(unique(unlist(strsplit(x = M_ID_a, split=","))))
  }
  M_c <- NULL
  if (!is.null(M_ID_c)) {
    M_c <- model_matrix_from_ID(M_ID_c)
    num_alleles_c <- length(unique(unlist(strsplit(x = M_ID_c, split=","))))
  }
  
  sim_data <- sim_mpp_qtl_partial(num_replicates = num_replicates,
                                  M_a = M_a,
                                  M_c = M_c,
                                  beta_a = beta_a,
                                  beta_c = beta_c,
                                  num_alleles_a = num_alleles_a, 
                                  num_alleles_c = num_alleles_c,
                                  qtl_a_effect_size = m_qtl_a_effect_size, 
                                  qtl_c_effect_size = y_qtl_c_effect_size, 
                                  m_b_effect_size = y_m_b_effect_size,
                                  sample_method = sample_method,
                                  num_founders = num_founders,
                                  m_strain_effect_size = m_strain_effect_size,
                                  y_strain_effect_size = y_strain_effect_size,
                                  m_noise_effect_size = m_noise_effect_size,
                                  y_noise_effect_size = y_noise_effect_size,
                                  impute = impute,
                                  locus_matrix = locus_matrix,
                                  return_means = return_means,
                                  num_sim = num_sim,
                                  sim_m_label = sim_m_label,
                                  sim_y_label = sim_y_label)
  
  return(list(m_data = sim_data$m_data,
              y_data = sim_data$y_data,
              locus_matrix = locus_matrix,
              properties=list(num_alleles_a = num_alleles_a,
                              num_alleles_c = num_alleles_c,
                              qtl_a_effect_size = original_effects$qtl_a_effect_size,
                              m_b_effect_size = original_effects$m_b_effect_size,
                              qtl_c_effect_size = original_effects$qtl_c_effect_size,
                              M_ID_a = M_ID_a,
                              M_ID_c = M_ID_c,
                              sample_method = sample_method,
                              num_replicates = num_replicates,
                              num_founders = num_founders,
                              strain_effect_size = original_effects$strain_effect_size,
                              impute = impute,
                              return_means = return_means)))
}

sim_mpp_qtl_partial <- function(locus_matrix,
                                num_replicates,
                                M_a = NULL,
                                M_c = NULL,
                                sample_method = c("uniform", "crp"),
                                qtl_a_effect_size, 
                                m_b_effect_size,
                                qtl_c_effect_size,
                                beta_a = NULL,
                                beta_c = NULL,
                                num_alleles_a = 8, 
                                num_alleles_c = 8, 
                                m_strain_effect_size,
                                y_strain_effect_size,
                                m_noise_effect_size,
                                y_noise_effect_size,
                                num_founders = 8,
                                num_sim,
                                impute = TRUE, 
                                return_means = TRUE,
                                sim_m_label = "sim_m",
                                sim_y_label = "sim_y",
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
  
  ## Strain effect 
  Z <- incidence_matrix(factor(strains, levels = strains))
  Z <- Z[rep(1:nrow(Z), each = num_replicates),]
  # For Y
  if (y_strain_effect_size != 0) {
    y_strain_effect <- rnorm(n = nrow(locus_matrix))
    y_strain_effect <- (y_strain_effect - mean(y_strain_effect))/sqrt(non_sample_var(y_strain_effect))
    y_strain_effect <- y_strain_effect * sqrt(y_strain_effect_size)
    
    y_strain_predictor <- Z %*% matrix(y_strain_effect, ncol = 1)
  }
  else { y_strain_predictor <- rep(0, nrow(Z)) }
  # For M
  if (m_strain_effect_size != 0) {
    m_strain_effect <- rnorm(n = nrow(locus_matrix))
    m_strain_effect <- (m_strain_effect - mean(m_strain_effect))/sqrt(non_sample_var(m_strain_effect))
    m_strain_effect <- m_strain_effect * sqrt(m_strain_effect_size)
    
    m_strain_predictor <- Z %*% matrix(m_strain_effect, ncol = 1)
  }
  else { m_strain_predictor <- rep(0, nrow(Z)) }
  
  #############
  ##
  ##   QTL
  ##
  #############
  ### Indirect effect (a) on M
  if (qtl_a_effect_size != 0) {
    qtl_a_effect <- sim_M_and_beta(num_alleles = num_alleles_a, 
                                   num_founders = num_founders, 
                                   M = M_a,
                                   sample_method = sample_method,
                                   beta = beta_a,
                                   ...)
  }
  else {
    qtl_a_effect <- list(M = diag(8),
                         beta = rep(0, 8))
  }
  M_a <- qtl_a_effect$M
  raw_beta_a <- qtl_a_effect$beta 
  if (qtl_a_effect_size != 0) {
    beta_a <- (raw_beta_a - mean(raw_beta_a))/sqrt(non_sample_var(raw_beta_a))
  }
  else { 
    beta_a <- raw_beta_a
  }
  ## Scaling to ZDMB, the observed variation in the population
  ZDMB_a <- Z %*% D %*% M_a %*% beta_a
  
  if (qtl_a_effect_size != 0) { # Case when more than one allele is observed
    qtl_a_predictor <- ZDMB_a * sqrt(qtl_a_effect_size/non_sample_var(ZDMB_a)[1])
  }
  else { # Case when only one allele is observed, should result in a null scan
    qtl_a_predictor <- ZDMB_a
  }
  
  ### Direct effect (c) on Y
  if (qtl_c_effect_size != 0) {
    qtl_c_effect <- sim_M_and_beta(num_alleles = num_alleles_c, 
                                   num_founders = num_founders, 
                                   M = M_c,
                                   sample_method = sample_method,
                                   beta = beta_c,
                                   ...)
  }
  else {
    qtl_c_effect <- list(M = diag(8),
                         beta = rep(0, 8))
  }
  M_c <- qtl_c_effect$M
  raw_beta_c <- qtl_c_effect$beta 
  if (qtl_c_effect_size != 0) {
    beta_c <- (raw_beta_c - mean(raw_beta_c))/sqrt(non_sample_var(raw_beta_c))
  }
  else { 
    beta_c <- raw_beta_c
  }
  # Scaling to ZDMB, the observed variation in the population
  ZDMB_c <- Z %*% D %*% M_c %*% beta_c
  
  if (qtl_c_effect_size != 0) { # Case when more than one allele is observed
    qtl_c_predictor <- ZDMB_c * sqrt(qtl_c_effect_size/non_sample_var(ZDMB_c)[1])
  }
  else { # Case when only one allele is observed, should result in a null scan
    qtl_c_predictor <- ZDMB_c
  }
  
  sim_m_data <- sim_y_data <- matrix(NA, nrow = nrow(Z), ncol = num_sim)
  for (i in 1:num_sim) {
    m_scaled_resid <- calc_scaled_residual(noise_effect_size = m_noise_effect_size,
                                           n = nrow(Z))
    y_scaled_resid <- calc_scaled_residual(noise_effect_size = y_noise_effect_size,
                                           n = nrow(Z))
    sim_m_data[,i] <- qtl_a_predictor + m_strain_predictor + m_scaled_resid
    
    m_predictor <- sim_m_data[,i] * sqrt(m_b_effect_size/non_sample_var(sim_m_data[,i])[1])
    
    sim_y_data[,i] <- m_predictor + qtl_c_predictor + y_strain_predictor + y_scaled_resid
  }
  colnames(sim_m_data) <- paste(sim_m_label, 1:ncol(sim_m_data), sep = "_")
  colnames(sim_y_data) <- paste(sim_y_label, 1:ncol(sim_y_data), sep = "_")
  rownames(sim_m_data) <- rownames(sim_y_data) <- paste(rep(strains, each = num_replicates), 1:num_replicates, sep = "_")
  
  if (num_replicates == 1) {
    rownames(sim_m_data) <- strains
    rownames(sim_y_data) <- strains
  } 
  else if (return_means) {
    sim_m_data <- data.frame(strain = rep(strains, each = num_replicates), sim_m_data) %>%
      tidyr::gather(key = "sim", value = "phenotype", -strain) %>%
      dplyr::group_by(strain, sim) %>%
      dplyr::summarize(phenotype = mean(phenotype)) %>%
      dplyr::ungroup() %>%
      tidyr::spread(key = "sim", value = "phenotype") %>%
      tibble::column_to_rownames("strain") %>%
      as.matrix
    sim_y_data <- data.frame(strain = rep(strains, each = num_replicates), sim_y_data) %>%
      tidyr::gather(key = "sim", value = "phenotype", -strain) %>%
      dplyr::group_by(strain, sim) %>%
      dplyr::summarize(phenotype = mean(phenotype)) %>%
      dplyr::ungroup() %>%
      tidyr::spread(key = "sim", value = "phenotype") %>%
      tibble::column_to_rownames("strain") %>%
      as.matrix
  }
  return(list(m_data = sim_m_data,
              y_data = sim_y_data,
              properties = list(qtl_a_effect_size = qtl_a_effect_size,
                                m_b_effect_size = m_b_effect_size,
                                qtl_c_effect_size = qtl_c_effect_size,
                                m_strain_effect_size = m_strain_effect_size,
                                y_strain_effect_size = y_strain_effect_size,
                                num_alleles_a = num_alleles_a,
                                num_alleles_c = num_alleles_c,
                                sample_method = sample_method,
                                num_replicates = num_replicates,
                                impute = impute,
                                return_means = return_means)))
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
#' @examples sim_target_from_mediator()
sim_target_from_mediator <- function(M,
                                     mediator_effect_size = 0.5,
                                     sim_label = "sim_y") {
  
  Y <- matrix(0, nrow = nrow(M), ncol = ncol(M))
  rownames(Y) <- rownames(M)
  colnames(Y) <- paste(sim_label, 1:ncol(M), sep = "_")
  for (i in 1:ncol(M)) {
    m_var <- non_sample_var(M[,i])
    y_var <- m_var * (1 - mediator_effect_size)
    e <- rnorm(n = nrow(M), mean = 0, sd = 1)
    if (mediator_effect_size != 0) {
      e_var <- non_sample_var(e)
      Y[,i] <- M[,i] * sqrt(mediator_effect_size/m_var)  + e * sqrt((1 - mediator_effect_size)/e_var)
    } else{
      Y[,i] <- M[,i] * sqrt(mediator_effect_size)  + e * sqrt(1 - mediator_effect_size)
    }
  }
  
  return(list(data = Y,
              mediator_effect_size = mediator_effect_size))
}

#' Simulate multiparental population (MPP) data from realized genomes
#'
#' This function takes qtl2 formatted genoprobs from an MPP population (e.g., CC or DO) and simulates 
#' outcomes driven by QTL.
#'
#' @param genoprobs A qtl2 formatted genoprobs object from qtl2. 
#' @param map A qtl2 formatted map object that corresponds to the genoprobs.
#' @param num_replicates DEFAULT: 1. The number of strain replicates to use.
#' @param num_sim DEFAULT: 1. The number of simulated outcomes to produce.
#' @param qtl_effect_size DEFAULT = 0.1. The proportion of variation in the simulated outcomes
#' explained by the QTL.
#' @export
#' @examples sim_mpp_data()
sim_mpp_data <- function(genoprobs, 
                         map,
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
  num_strains <- dim(genoprobs[[1]])[1]
  
  original_effects <- list(qtl_effect_size = qtl_effect_size,
                           strain_effect_size = strain_effect_size)
  
  ############ Scaling effects
  noise_effect_size <- (1 - qtl_effect_size - strain_effect_size)
  reduced_noise_effect_size <- noise_effect_size/num_replicates
  scale_factor <- 1/sum(c(qtl_effect_size, strain_effect_size, reduced_noise_effect_size))
  qtl_effect_size <- scale_factor * qtl_effect_size
  strain_effect_size <- scale_factor * strain_effect_size
  noise_effect_size <- scale_factor * noise_effect_size
  
  ## Number of individuals
  num_ind <- ifelse(return_means, num_strains, num_strains * num_replicates)
  
  ## Setting number of alleles to match pre-specified QTL effect - convenience
  if (!is.null(beta)) {
    num_alleles <- length(beta)
  }
  
  ## Sampling loci
  loci <- unlist(lapply(genoprobs, function(x) dimnames(x)[[3]]))
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
    this_sim <- sim_mpp_qtl(num_replicates = num_replicates,
                            M = M,
                            beta = beta,
                            sample_method = sample_method,
                            num_alleles = num_alleles, 
                            num_founders = num_founders,
                            qtl_effect_size = qtl_effect_size, 
                            strain_effect_size = strain_effect_size,
                            noise_effect_size = noise_effect_size,
                            impute = impute,
                            locus_matrix = qtl2::pull_genoprobpos(genoprobs = genoprobs, marker = locus[locus_index[i]]),
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
  
  map_df <- qtl2convert::map_list_to_df(map)

  return(list(data = sim_matrix,
              locus = as.character(locus),
              locus_pos = as.numeric(map_df[locus, "pos"]),
              locus_chr = as.character(map_df[locus, "chr"]),
              properties = list(num_alleles = num_alleles,
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

#' Simulate multiparental population (MPP) data from a specified locus design matrix
#'
#' This function takes a locus design matrix, e.g., founder haplotypes dosages from CC mice,
#' to simulate outcomes driven by its QTL.
#'
#' @param locus_matrix A design matrix of founder haplotypes at QTL for MPP founder strains. The number of 
#' rows should correspond to the number of strains. The number of columns is eight for the CC/DO. 
#' @param num_replicates DEFAULT: 1. The number of strain replicates to use.
#' @param num_sim DEFAULT: 1. The number of simulated outcomes to produce.
#' @param qtl_effect_size DEFAULT = 0.1. The proportion of variation in the simulated outcomes
#' explained by the QTL.
#' @export
#' @examples sim_cc_single_locus()
sim_mpp_single_locus <- function(locus_matrix, 
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
  scale_factor <- 1/sum(c(qtl_effect_size, strain_effect_size, reduced_noise_effect_size))
  qtl_effect_size <- scale_factor * qtl_effect_size
  strain_effect_size <- scale_factor * strain_effect_size
  noise_effect_size <- scale_factor * noise_effect_size
  
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
  
  sim_data <- sim_mpp_qtl(num_replicates = num_replicates,
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
                          num_sim = num_sim,
                          sim_label = sim_label)
  
  return(list(data = sim_data$data,
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

sim_mpp_qtl <- function(locus_matrix,
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
    qtl_effect <- sim_M_and_beta(num_alleles = num_alleles, 
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
  
  ## Scaling to ZDMB, the observed variation in the population
  ZDMB <- Z %*% D %*% M %*% beta
  
  if (qtl_effect_size != 0) { # Case when more than one allele is observed
    qtl_predictor <- ZDMB * sqrt(qtl_effect_size/non_sample_var(ZDMB)[1])
  }
  else { # Case when only one allele is observed, should result in a null scan
    qtl_predictor <- ZDMB
  }
  
  sim_data <- matrix(NA, nrow = nrow(Z), ncol = num_sim)
  for (i in 1:num_sim) {
    scaled_resid <- calc_scaled_residual(noise_effect_size = noise_effect_size,
                                         n = nrow(Z))
    sim_data[,i] <- qtl_predictor + strain_predictor + scaled_resid
  }

  colnames(sim_data) <- paste(sim_label, 1:ncol(sim_data), sep = "_")
  rownames(sim_data) <- paste(rep(strains, each = num_replicates), 1:num_replicates, sep = "_")
  
  if (num_replicates == 1) {
    rownames(sim_data) <- strains
  } 
  else if (return_means) {
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
sim_M_and_beta <- function(num_alleles = 8, 
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

## Draws and scales residuals in single function
calc_scaled_residual <- function(noise_effect_size, n) {
  
  residual <- rnorm(n = n)
  residual <- (residual - mean(residual))/sqrt(non_sample_var(residual))
  residual <- residual*sqrt(noise_effect_size)
  residual
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



