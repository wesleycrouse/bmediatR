#' Sum-to-zero contrast matrix
#'
#' This function takes the dimensions of a design matrix to produce a sum-to-zero contrast matrix, which is post-multiplied
#' by the design matrix, reducing its dimension by one, resulting in a parameterization of the model in terms of deviation
#' from mean (of a balanced data set). 
#'
#' @param k Number of columns of the design matrix.
#' @export
#' @examples sumtozero_contrast()
sumtozero_contrast <- function(k) {
  u <- 1/((k-1)^(0.5))
  v <- (-1+(k^(0.5)))*((k-1)^(-1.5))
  w <- (k-2)*v+u
  C <- matrix(-v, k-1, k-1)
  diag(C) <- rep(w, k-1)
  rbind(C, rep(-u, k-1))
}

# Adapted from mvtnorm::dmvt
dmvt_chol <- function(x, sigma_chol, df) {
  if (is.vector(x)) { x <- matrix(x, nrow=length(x)) }
  n <- nrow(x)
  dec <- sigma_chol
  R_x_m <- backsolve(dec, x, transpose = TRUE)
  rss <- colSums(R_x_m^2)
  lgamma(0.5*(df+n)) - (lgamma(0.5*df) + sum(log(diag(dec))) + 0.5*n*log(pi*df)) - 0.5*(df+n)*log1p(rss/df)
}

# Adapted from qtl2::batch_cols
batch_cols <- function(mat) {
  mat <- !is.finite(mat)
  n <- nrow(mat)
  all_true <- rep(TRUE, n)
  result <- NULL
  n_na <- colSums(mat)
  no_na <- (n_na == 0)
  if (any(no_na)){ result <- list(list(cols=which(no_na), omit=numeric(0))) }
  one_na <- (n_na == 1)
  if (any(one_na)) {
    wh <- apply(mat[,one_na, drop = FALSE], 2, which)
    spl <- split(which(one_na), wh)
    part2 <- lapply(seq_along(spl), function(i){ list(cols=as.numeric(spl[[i]]), omit=as.numeric(names(spl)[i])) })
    if (is.null(result)) {
      result <- part2
    } else {
      result <- c(result, part2)
    }
  }
  other_cols <- !(no_na | one_na)
  if (any(other_cols)) {
    other_cols <- seq_len(ncol(mat))[other_cols]
    pat <- apply(mat[,other_cols, drop = FALSE], 2, function(a) paste(which(a), collapse = ":"))
    u <- unique(pat)
    part3 <- lapply(u, function(a) list(cols=other_cols[pat==a], omit=as.numeric(strsplit(a, ":")[[1]])))
  } else {
    part3 <- NULL
  }
  if (is.null(result)) {
    part3
  } else {
    c(result, part3)
  }
}

#' Convert preset options to log prior case pobabilities function 
#'
#' This function takes the log prior case probabilities, and if a preset is provided, converts it into the formal log prior case
#' probability.
#'
#' @param ln_prior_c Log prior case probabilities. If posterior_summary() is being used for a non-default posterior odds
#' summary, the log prior case probabilities used with bmediatR() are stored in its output.
#' @export
#' @examples return_ln_prior_c_from_presets()
return_ln_prior_c_from_presets <- function(ln_prior_c) {
  
  if (ln_prior_c[1] == "complete"){
    ln_prior_c <- c(rep(0,8), rep(-Inf,4))
  } else if (ln_prior_c[1] == "partial"){
    ln_prior_c <- c(rep(-Inf,4), rep(0,4), rep(-Inf,4))
  } else if (ln_prior_c[1] == "reactive"){
    ln_prior_c <- rep(0,12)
  } else {
    ln_prior_c <- ln_prior_c
  }
  ln_prior_c
}

#' Summaries of posterior model probabilities function 
#'
#' This function takes the log posterior probability of the data (posterior likelihood) for the various cases, the log prior case probabilities, and
#' returns log posterior odds.
#'
#' @param ln_prob_data Log posterior likelihoods under the various models, returned by bmediatR().
#' @param ln_prior_c Log prior case probabilities. If posterior_summary() is being used for a non-default posterior odds
#' summary, the log prior case probabilities used with bmediatR() are stored in its output.
#' @param c_numerator The index of cases to be summed in the numerator of the posterior odds. Cases, their order, and likelihoods
#' are provided in model_info().
#' @export
#' @examples posterior_summary()
posterior_summary <- function(ln_prob_data, 
                              ln_prior_c, 
                              c_numerator){
  
  #function to compute log odds from log probabilities
  ln_odds <- function(ln_p, numerator){
    ln_odds_numerator <- apply(ln_p[,numerator, drop=F], 1, matrixStats::logSumExp)
    ln_odds_denominator <- apply(ln_p[,-numerator, drop=F], 1, matrixStats::logSumExp)
    ln_odds <- ln_odds_numerator -ln_odds_denominator
  }
  
  #ensure c_numerator is a list
  if (!is.list(c_numerator)){
    c_numerator <- list(c_numerator)
  }
  
  #presets for ln_prior_c; 
  ln_prior_c <- return_ln_prior_c_from_presets(ln_prior_c = ln_prior_c)
  
  #ensure ln_prior_c sum to 1 on probability scale and that it is a matrix
  if (is.matrix(ln_prior_c)){
    ln_prior_c <- t(apply(ln_prior_c, 1, function(x){x - matrixStats::logSumExp(x)}))
  } else {
    ln_prior_c <- ln_prior_c - matrixStats::logSumExp(ln_prior_c)
    ln_prior_c <- matrix(ln_prior_c, nrow(ln_prob_data), length(ln_prior_c), byrow=T)
  }
  
  #compute posterior probabilities for all cases
  #cases encoded by presence (1) or absence (0) of 'X->y, X->m, m->y' edges on the DAG
  #(*) denotes reverse causation 'm<-y'
  #c1:  '0,0,0' / H1 and H5
  #c2:  '0,0,1' / H2 and H5
  #c3:  '0,1,0' / H1 and H6
  #c4:  '0,1,1' / H2 and H6 - complete mediation
  #c5:  '1,0,0' / H3 and H5 
  #c6:  '1,0,1' / H4 and H5
  #c7:  '1,1,0' / H3 and H6 - colocalization
  #c8:  '1,1,1' / H4 and H6 - partial mediation
  #c9:  '0,0,*' / H1 and H7
  #c10: '0,1,*' / H1 and H8
  #c11: '1,0,*' / H3 and H7
  #c12: '1,1,*' / H3 and H8
  ln_post_c <- cbind(ln_prob_data[,1] + ln_prob_data[,5] + ln_prior_c[,1],
                     ln_prob_data[,2] + ln_prob_data[,5] + ln_prior_c[,2],
                     ln_prob_data[,1] + ln_prob_data[,6] + ln_prior_c[,3],
                     ln_prob_data[,2] + ln_prob_data[,6] + ln_prior_c[,4],
                     ln_prob_data[,3] + ln_prob_data[,5] + ln_prior_c[,5],
                     ln_prob_data[,4] + ln_prob_data[,5] + ln_prior_c[,6],
                     ln_prob_data[,3] + ln_prob_data[,6] + ln_prior_c[,7],
                     ln_prob_data[,4] + ln_prob_data[,6] + ln_prior_c[,8],
                     ln_prob_data[,1] + ln_prob_data[,7] + ln_prior_c[,9],
                     ln_prob_data[,1] + ln_prob_data[,8] + ln_prior_c[,10],
                     ln_prob_data[,3] + ln_prob_data[,7] + ln_prior_c[,11],
                     ln_prob_data[,3] + ln_prob_data[,8] + ln_prior_c[,12])
  
  ln_ml <- apply(ln_post_c, 1, matrixStats::logSumExp)
  ln_post_c <- ln_post_c - ln_ml
  
  colnames(ln_post_c) <- c("0,0,0",
                           "0,0,1",
                           "0,1,0",
                           "0,1,1",
                           "1,0,0",
                           "1,0,1",
                           "1,1,0",
                           "1,1,1",
                           "0,0,*",
                           "0,1,*",
                           "1,0,*",
                           "1,1,*")
  rownames(ln_post_c) <- rownames(ln_prob_data)
  
  #compute prior odds for each combination of cases
  ln_prior_odds <- sapply(c_numerator, ln_odds, ln_p=ln_prior_c)
  ln_prior_odds <- matrix(ln_prior_odds, ncol=length(c_numerator))
  rownames(ln_prior_odds) <- rownames(ln_post_c)

  #compute posterior odds for each combination of cases
  ln_post_odds <- sapply(c_numerator, ln_odds, ln_p=ln_post_c)
  ln_post_odds <- matrix(ln_post_odds, ncol=length(c_numerator))
  rownames(ln_post_odds) <- rownames(ln_post_c)
  
  if (is.null(c_numerator)) {
    colnames(ln_post_odds) <- colnames(ln_prior_odds) <- c_numerator
  } else {
    colnames(ln_post_odds) <- colnames(ln_prior_odds) <- names(c_numerator)
  }
  
  #return results
  list(ln_post_c=ln_post_c, ln_post_odds=ln_post_odds, ln_prior_odds=ln_prior_odds, ln_ml=ln_ml)
}

#' Column indeces for commonly used posterior odds
#'
#' This helper function returns the columns of the log posterior case probabilities to be summed for
#' commonly desired log posterior odds summaries.
#'
#' @param odds_type The desired posterior odds.
#' @export
#' @examples return_preset_odds_index()
return_preset_odds_index <- function(odds_type = c("mediation", 
                                                   "partial", 
                                                   "complete", 
                                                   "colocal", 
                                                   "mediation_or_colocal",
                                                   "y_depends_x", 
                                                   "reactive")) {
  
  presets <- list("mediation" = c(4, 8),
                  "partial" = 8,
                  "complete" = 4,
                  "colocal" = 7,
                  "mediation_or_colocal" = c(4, 7, 8),
                  "y_depends_x" = c(4:8, 11, 12),
                  "reactive" = 9:12)
  
  index_list <- presets[odds_type]
  index_list
}

## Function to process data and optionally align them
process_data <- function(y, M, X, 
                         Z = NULL, Z_y = NULL, Z_M = NULL,
                         w = NULL, w_y = NULL, w_M = NULL,
                         align_data = TRUE,
                         verbose = TRUE) {
  
  # Ensure y is a vector
  if (is.matrix(y)) { y <- y[,1] }

  # Ensure X, M, Z, Z_y, and Z_M are matrices
  X <- as.matrix(X)
  M <- as.matrix(M)
  if (!is.null(Z)) { Z <- as.matrix(Z) }
  if (!is.null(Z_y)) { Z_y <- as.matrix(Z_y) }
  if (!is.null(Z_M)) { Z_M <- as.matrix(Z_M) }
  
  # Process covariate matrices
  if (is.null(Z_y)) { Z_y <- matrix(NA, length(y), 0); rownames(Z_y) <- names(y) }
  if (is.null(Z_M)) { Z_M <- matrix(NA, nrow(M), 0); rownames(Z_M) <- rownames(M) }
  
  if (!is.null(Z)) {
    if (align_data) {
      Z_y <- cbind(Z, Z_y[rownames(Z),])
      Z_M <- cbind(Z, Z_M[rownames(Z),])
    } else {
      Z_y <- cbind(Z, Z_y)
      Z_M <- cbind(Z, Z_M)
    }
  }
  
  # Process weight vectors
  if (is.null(w)) {
    if (is.null(w_y)) { w_y <- rep(1, length(y)); names(w_y) <- names(y) }
    if (is.null(w_M)) { w_M <- rep(1, nrow(M)); names(w_M) <- rownames(M) }
  } else {
    w_y <- w_M <- w
  }
  
  if (align_data) {
    # M and Z_M can have NAs
    overlapping_samples <- Reduce(f = intersect, x = list(names(y), 
                                                          rownames(X),
                                                          rownames(Z_y), 
                                                          names(w_y)))
    
    if (length(overlapping_samples) == 0 | !any(overlapping_samples %in% unique(c(rownames(M), rownames(Z_M), names(w_M))))) {
      stop("No samples overlap. Check rownames of M, X, Z (or Z_y and Z_M) and names of y and w (or w_y and w_M).", call. = FALSE)
    } else if (verbose) {
      writeLines(text = c("Number of overlapping samples:", length(overlapping_samples)))
    }
    # Ordering
    y <- y[overlapping_samples]
    M <- M[overlapping_samples,, drop = FALSE]
    X <- X[overlapping_samples,, drop = FALSE]
    Z_y <- Z_y[overlapping_samples,, drop = FALSE]
    Z_M <- Z_M[overlapping_samples,, drop = FALSE]
    w_y <- w_y[overlapping_samples]
    w_M <- w_M[overlapping_samples]
  }
  
  # Drop observations with missing y or X and update n
  complete_y <- !is.na(y)
  complete_X <- !apply(is.na(X), 1, any)
  
  y <- y[complete_y & complete_X]
  M <- M[complete_y & complete_X,, drop = FALSE]
  X <- X[complete_y & complete_X,, drop = FALSE]
  Z_y <- Z_y[complete_y & complete_X,, drop = FALSE]
  Z_M <- Z_M[complete_y & complete_X,, drop = FALSE]
  w_y <- w_y[complete_y & complete_X]
  w_M <- w_M[complete_y & complete_X]
  
  # Drop columns of Z_y and Z_M that are invariant
  Z_y_drop <- which(apply(Z_y, 2, function(x) var(x)) == 0)
  Z_M_drop <- which(apply(Z_M, 2, function(x) var(x)) == 0)
  if (length(Z_y_drop) > 0) {
    writeLines(paste("Dropping invariants columns from Z_y:", colnames(Z_y)[Z_y_drop]))
    Z_y <- Z_y[,-Z_y_drop, drop = FALSE]
  }
  if (length(Z_M_drop) > 0) {
    writeLines(paste("Dropping invariants columns from Z_M:", colnames(Z_M)[Z_M_drop]))
    Z_M <- Z_M[,-Z_M_drop, drop = FALSE]
  }

  # Return aligned data
  list(y = y,
       M = M,
       X = X,
       Z_y = Z_y, Z_M = Z_M,
       w_y = w_y, w_M = w_M)
}

#' Bayesian model selection for mediation analysis function 
#'
#' This function takes an outcome (y), candidate mediators (M), and a driver as a design matrix (X) to perform a 
#' Bayesian model selection analysis for mediation.
#' 
#' @param y Vector or single column matrix of an outcome variable. Single outcome variable expected. 
#' Names or rownames must match across M, X, Z, Z_y, Z_M, w, w_y, and w_M (if provided) when align_data = TRUE. If align_data = FALSE,
#' dimensions and order must match across inputs.
#' @param M Vector or matrix of mediator variables. Multiple mediator variables are supported. 
#' Names or rownames must match across y, X, Z, Z_y, Z_M, w, w_y, and w_M (if provided) when align_data = TRUE. If align_data = FALSE,
#' dimensions and order must match across inputs.
#' @param X Design matrix of the driver. Names or rownames must match across y, M, Z, Z_y, Z_M, w, w_y, and w_M (if provided) when align_data = TRUE. 
#' If align_data = FALSE, dimensions and order must match across inputs. One common application is for X to represent genetic information at a QTL, 
#' either as founder strain haplotypes or variant genotypes, though X is generalizable to other types of variables.
#' @param Z DEFAULT: NULL. Design matrix of covariates that influence the outcome and mediator variables. 
#' Names or rownames must match to those of y, M, X, w, w_y, and w_M (if provided) when align_data = TRUE. If align_data=FALSE,
#' dimensions and order must match across inputs. If Z is provided, it is added interally to Z_y and Z_M.
#' @param Z_y DEFAULT: NULL. Design matrix of covariates that influence the outcome variable. 
#' Names or rownames must match to those of y, M, X, Z_M, w, w_y, and w_M (if provided) when align_data = TRUE. If align_data = FALSE,
#' dimensions and order must match across inputs.
#' @param Z_M DEFAULT: NULL. Design matrix of covariates that influence the mediator variables. 
#' Names or rownames must match across y, M, X, Z_y, w, w_y, and w_M (if provided) when align_data = TRUE. If align_data = FALSE,
#' dimensions and order must match across inputs.
#' @param w DEFAULT: NULL. Vector or single column matrix of weights for individuals in analysis that applies to both 
#' y and M. Names must match across y, M, X, Z, Z_y, and Z_M (if provided) when align_data = TRUE. If align_data = FALSE,
#' dimensions and order must match across inputs. A common use would be for an analysis of strain means, where w 
#' is a vector of the number of individuals per strain. If no w, w_y, or w_M is given, observations are equally weighted as 1s for y and M. 
#' If w is provided, it supercedes w_y and w_M.
#' @param w_y DEFAULT: NULL. Vector or single column matrix of weights for individuals in analysis, specific to the measurement
#' of y. Names must match across y, M, X, Z, Z_y, Z_M, and w_M (if provided) when align_data = TRUE. If align_data = FALSE,
#' dimensions and order must match across inputs. A common use would be for an analysis of strain means, where y and M
#' are summarized from a different number of individuals per strain. w_y is a vector of the number of individuals per strain used to
#' measure y. If no w_y (or w) is given, observations are equally weighted as 1s for y.
#' @param w_M DEFAULT: NULL. Vector or single column matrix of weights for individuals in analysis, specific to the measurement 
#' of M. Names must match across y, M, X, Z, Z_y, Z_M, and w_y (if provided) when align_data = TRUE. If align_data = FALSE,
#' dimensions and order must match across inputs. A common use would be for an analysis of strain means, where y and M
#' are summarized from a different number of individuals per strain. w_M is a vector of the number of individuals per strain use to 
#' measure M. If no w_M (or w) is given, observations are equally weighted as 1s for M.
#' @param tau_sq_mu DEFAULT: 1000. Variance component for the intercept. The DEFAULT represents a diffuse prior, analagous to 
#' a fixed effect term.
#' @param tau_sq_Z DEFAULT: 1000. Variance component for the covariates encoded in Z. The DEFAULT represents a diffuse prior, analagous 
#' to fixed effect terms.
#' @param phi_sq DEFAULT: c(1, 1, 1). Each element of (a, b, c) represents one of the relationships being evaluated for mediation, 
#' specifically the ratio of signal to noise. a is the effect of X on M, b is the effect of M on y, and c is the effect of X on y. 
#' The DEFAULT represents relationships that explain 50% of the variation in the outcome variable.
#' @param ln_prior_c DEFAULT: "complete". The prior log case probabilities. See model_info() for description of likelihoods and their
#' combinations into cases. Simplified pre-set options are available, including "complete", "partial", and "reactive".
#' @param align_data DEFAULT: TRUE. If TRUE, expect vector and matrix inputes to have names and rownames, respectively. The overlapping data
#' will then be aligned, allowing the user to not have to reduce data to overlapping samples and order them.
#' @export
#' @examples bmediatR()
bmediatR <- function(y, M, X, 
                     Z = NULL, Z_y = NULL, Z_M = NULL,
                     w = NULL, w_y = NULL, w_M = NULL,
                     kappa = 0.001, lambda = 0.001,
                     tau_sq_mu = 1000, tau_sq_Z = 1000,
                     phi_sq = c(1, 1, 1),
                     ln_prior_c = "complete",
                     options_X = list(sum_to_zero = TRUE, center = FALSE, scale = FALSE),
                     align_data = TRUE,
                     verbose = TRUE) {
  
  #presets for ln_prior_c; 
  ln_prior_c <- return_ln_prior_c_from_presets(ln_prior_c = ln_prior_c)
  
  #ensure ln_prior_c sum to 1 on probability scale
  if (is.matrix(ln_prior_c)){
    ln_prior_c <- t(apply(ln_prior_c, 1, function(x){x - matrixStats::logSumExp(x)}))
  } else {
    ln_prior_c <- ln_prior_c - matrixStats::logSumExp(ln_prior_c)
  }

  #optionally align data
  processed_data <- process_data(y = y, M = M, X = X,
                                 Z_y = Z_y, Z_M = Z_M,
                                 w_y = w_y, w_M = w_M, 
                                 align_data = align_data,
                                 verbose = verbose)
  y <- processed_data$y
  M <- processed_data$M
  X <- processed_data$X
  Z_y <- processed_data$Z_y
  Z_M <- processed_data$Z_M
  w_y <- processed_data$w_y
  w_M <- processed_data$w_M
  
  #dimenion of y
  n <- length(y)
  
  #dimension of Z's
  p_y <- ncol(Z_y)
  p_M <- ncol(Z_M)
  
  #scale y, M, and Z
  y <- c(scale(y))
  M <- apply(M, 2, scale)
  if (p_y > 0) { Z_y <- apply(Z_y, 2, scale) }
  if (p_M > 0) { Z_M <- apply(Z_M, 2, scale) }
  
  #optionally use sum-to-zero contrast for X
  #recommended when X is a matrix of factors, with a column for every factor level
  if (options_X$sum_to_zero == TRUE) {
    C <- sumtozero_contrast(ncol(X))
    X <- X%*%C
  }
  
  #optionally center and scale X
  X <- apply(X, 2, scale, center = options_X$center, scale = options_X$scale)
  
  #dimension of X
  d <- ncol(X)
  
  #column design matrix for mu
  ones <- matrix(1, n)
  
  #begin Bayesian calculations
  if (verbose) { print("Initializing", quote = FALSE) }
  
  #reformat priors
  kappa = rep(kappa, 8)
  lambda = rep(lambda, 8)
  tau_sq_mu = rep(tau_sq_mu, 8)
  tau_sq_Z = rep(tau_sq_Z, 8)
  phi_sq_X = c(NA,NA,phi_sq[3],phi_sq[3],NA,phi_sq[1],NA,phi_sq[1])
  phi_sq_m = c(NA,phi_sq[2],NA,phi_sq[2],NA,NA,NA,NA)
  phi_sq_y = c(NA,NA,NA,NA,NA,NA,phi_sq[2],phi_sq[2])
  
  #identify likelihoods that are not supported by the prior
  #will not compute cholesky or likelihood for these
  calc_ln_prob_data <- rep(NA, 8) 
  calc_ln_prob_data[1] <- any(!is.infinite(ln_prior_c[c(1,3,9,10)]))
  calc_ln_prob_data[2] <- any(!is.infinite(ln_prior_c[c(2,4)]))
  calc_ln_prob_data[3] <- any(!is.infinite(ln_prior_c[c(5,7,11,12)]))
  calc_ln_prob_data[4] <- any(!is.infinite(ln_prior_c[c(6,8)]))
  calc_ln_prob_data[5] <- any(!is.infinite(ln_prior_c[c(1,2,5,6)]))
  calc_ln_prob_data[6] <- any(!is.infinite(ln_prior_c[c(3,4,7,8)]))
  calc_ln_prob_data[7] <- any(!is.infinite(ln_prior_c[c(9,11)]))
  calc_ln_prob_data[8] <- any(!is.infinite(ln_prior_c[c(10,12)]))
  
  #likelihood models for all hypothesis
  #hypotheses encoded by presence (1) or absence (0) of 'X->y, X->m, m->y' edges on the DAG
  #(*) denotes reverse causation 'm<-y'
  #H1: '0,-,0' / y does not depend on X or m
  #H2: '0,-,1' / y depends on m but not X
  #H3: '1,-,0' / y depends on X but not m
  #H4: '1,-,1' / y depends on X and m
  #H5: '-,0,-' / m does not depend on X
  #H6: '-,1,-' / m depends on X
  #H7: '-,0,*' / m depends on y but not X
  #H8: '-,1,*' / m depends on X and y
  #all include covariates Z
  
  #design matrices for H1,H3,H5-H8 complete cases (do not depend on m)
  X1 <- cbind(ones, Z_y)
  X3 <- cbind(ones, X, Z_y)
  X5 <- cbind(ones, Z_M)
  X6 <- cbind(ones, X, Z_M)
  X7 <- cbind(X5, y)
  X8 <- cbind(X6, y)
  
  #check if all scale hyperparameters are identical for H1 and H5
  #implies sigma1 and sigma5 identical, used to reduce computations
  sigma5_equal_sigma1 <- all(lambda[1]==lambda[5],
                             tau_sq_mu[1] == tau_sq_mu[5], 
                             tau_sq_Z[1] == tau_sq_Z[5],
                             identical(Z_y, Z_M),
                             identical(w_y, w_M))
  
  #check if all scale hyperparameters are identical for H3 and H6
  #implies sigma3 and sigma6 identical, used to reduce computations
  sigma6_equal_sigma3 <- all(lambda[3]==lambda[6],
                             tau_sq_mu[3] == tau_sq_mu[6], 
                             tau_sq_Z[3] == tau_sq_Z[6],
                             identical(Z_y, Z_M),
                             identical(w_y, w_M))
  
  #prior variance matrices (diagonal) for H1-H8 
  v1 <- c(tau_sq_mu[1], rep(tau_sq_Z[1], p_y))
  v2 <- c(tau_sq_mu[2], rep(tau_sq_Z[2], p_y), phi_sq_m[2])
  v3 <- c(tau_sq_mu[3], rep(phi_sq_X[3], d), rep(tau_sq_Z[3], p_y))
  v4 <- c(tau_sq_mu[4], rep(phi_sq_X[4], d), rep(tau_sq_Z[4], p_y), phi_sq_m[4])
  v7 <- c(tau_sq_mu[7], rep(tau_sq_Z[7], p_M), phi_sq_y[7])
  v8 <- c(tau_sq_mu[8], rep(phi_sq_X[8], d), rep(tau_sq_Z[8], p_M), phi_sq_y[8])
  
  if (!sigma5_equal_sigma1 | !calc_ln_prob_data[1]) {
    v5 <- c(tau_sq_mu[5], rep(tau_sq_Z[5], p_M))
  }
  
  if (!sigma6_equal_sigma3 | !calc_ln_prob_data[3]) {
    v6 <- c(tau_sq_mu[6], rep(phi_sq_X[6], d), rep(tau_sq_Z[6], p_M))
  }
  
  #scale matrices for H1,H3,H5-H8 complete cases (do not depend on m)
  sigma1 <- crossprod(sqrt(lambda[1]*v1)*t(X1))
  sigma3 <- crossprod(sqrt(lambda[3]*v3)*t(X3))
  sigma7 <- crossprod(sqrt(lambda[7]*v7)*t(X7))
  sigma8 <- crossprod(sqrt(lambda[8]*v8)*t(X8))
  
  diag(sigma1) <- diag(sigma1) + lambda[1]/w_y
  diag(sigma3) <- diag(sigma3) + lambda[3]/w_y
  diag(sigma7) <- diag(sigma7) + lambda[7]/w_M
  diag(sigma8) <- diag(sigma8) + lambda[8]/w_M
  
  if (!sigma5_equal_sigma1 | !calc_ln_prob_data[1]) {
    sigma5 <- crossprod(sqrt(lambda[5]*v5)*t(X5))
    diag(sigma5) <- diag(sigma5) + lambda[5]/w_M
  }
  
  if (!sigma6_equal_sigma3 | !calc_ln_prob_data[3]) {
    sigma6 <- crossprod(sqrt(lambda[6]*v6)*t(X6))
    diag(sigma6) <- diag(sigma6) + lambda[6]/w_M
  }
  
  #object to store likelihoods
  ln_prob_data=matrix(-Inf, ncol(M), 8)
  rownames(ln_prob_data) <- colnames(M)
  colnames(ln_prob_data) <- c("0,-,0", 
                              "0,-,1", 
                              "1,-,0", 
                              "1,-,1",
                              "-,0,-", 
                              "-,1,-",
                              '-,0,*',
                              '-,1,*')
  
  #identify batches of M that have the same pattern of missing values
  missing_m <- bmediatR:::batch_cols(M)
  
  #iterate over batches of M with same pattern of missing values
  if(verbose){print("Iterating", quote=F)}
  counter <- 0
  
  for (b in 1:length(missing_m)){
    #subset to non-missing observations
    index <- rep(T, length(y))
    index[missing_m[[b]]$omit] <- F
    
    if (any(index)){
      y_subset <- y[index]
      w_y_subset <- w_y[index]
      w_M_subset <- w_M[index]
      
      #cholesky matrices for H1,H3,H5-H8 non-missing observations (do not depend on m)
      if (calc_ln_prob_data[1]) { sigma1_chol_subset <- chol(sigma1[index,index]) }
      if (calc_ln_prob_data[3]) { sigma3_chol_subset <- chol(sigma3[index,index]) }
      if (calc_ln_prob_data[7]) { sigma7_chol_subset <- chol(sigma7[index,index]) }
      if (calc_ln_prob_data[8]) { sigma8_chol_subset <- chol(sigma8[index,index]) }
      
      if (sigma5_equal_sigma1 & calc_ln_prob_data[1]) {
        sigma5_chol_subset <- sigma1_chol_subset
      } else if (calc_ln_prob_data[5]) {
        sigma5_chol_subset <- chol(sigma5[index,index])
      }
      
      if (sigma6_equal_sigma3 & calc_ln_prob_data[3]) {
        sigma6_chol_subset <- sigma3_chol_subset
      } else if (calc_ln_prob_data[6]) {
        sigma6_chol_subset <- chol(sigma6[index,index])
      }
      
      #compute H1 and H3 outside of the mediator loop (invariant)
      if (calc_ln_prob_data[1]) { ln_prob_data1 <- bmediatR:::dmvt_chol(y_subset, sigma_chol=sigma1_chol_subset, df = kappa[1]) }
      if (calc_ln_prob_data[3]) { ln_prob_data3 <- bmediatR:::dmvt_chol(y_subset, sigma_chol=sigma3_chol_subset, df = kappa[3]) }
      
      #iterate over mediators
      for (i in missing_m[[b]]$cols) {
        counter <- counter + 1
        if (counter%%1000==0 & verbose) { print(paste(counter, "of", ncol(M)), quote=F) }
        
        #set current mediator non-missing observations
        m_subset <- M[index,i]
        
        #design matrix for H2 and H4 non-missing observations
        X2_subset <- cbind(X1[index,,drop=F], m_subset)
        X4_subset <- cbind(X3[index,,drop=F], m_subset)
        
        #scale and cholesky matrices for H2 and H4 non-missing observations
        sigma2_subset <- crossprod(sqrt(lambda[2]*v2)*t(X2_subset))
        sigma4_subset <- crossprod(sqrt(lambda[4]*v4)*t(X4_subset))
        
        diag(sigma2_subset) <- diag(sigma2_subset) + lambda[2]/w_y_subset
        diag(sigma4_subset) <- diag(sigma4_subset) + lambda[4]/w_y_subset
        
        if (calc_ln_prob_data[2]) { sigma2_chol_subset <- chol(sigma2_subset) }
        if (calc_ln_prob_data[4]) { sigma4_chol_subset <- chol(sigma4_subset) }
        
        #compute likelihoods for H1-H8
        if (calc_ln_prob_data[1]) { ln_prob_data[i,1] <- ln_prob_data1 }
        if (calc_ln_prob_data[3]) { ln_prob_data[i,3] <- ln_prob_data3 }
        
        if (calc_ln_prob_data[2]) { ln_prob_data[i,2] <- bmediatR:::dmvt_chol(y_subset, sigma_chol=sigma2_chol_subset, df = kappa[2]) }
        if (calc_ln_prob_data[4]) { ln_prob_data[i,4] <- bmediatR:::dmvt_chol(y_subset, sigma_chol=sigma4_chol_subset, df = kappa[4]) }
        if (calc_ln_prob_data[5]) { ln_prob_data[i,5] <- bmediatR:::dmvt_chol(m_subset, sigma_chol=sigma5_chol_subset, df = kappa[5]) }
        if (calc_ln_prob_data[6]) { ln_prob_data[i,6] <- bmediatR:::dmvt_chol(m_subset, sigma_chol=sigma6_chol_subset, df = kappa[6]) }
        if (calc_ln_prob_data[7]) { ln_prob_data[i,7] <- bmediatR:::dmvt_chol(m_subset, sigma_chol=sigma7_chol_subset, df = kappa[7]) }
        if (calc_ln_prob_data[8]) { ln_prob_data[i,8] <- bmediatR:::dmvt_chol(m_subset, sigma_chol=sigma8_chol_subset, df = kappa[8]) }
      }
    }
  }
  
  #compute posterior probabilities for all cases
  #compute posterior odds for specified combinations of cases
  #cases encoded by presence (1) or absence (0) of 'X->y, X->m, m->y' edges on the DAG
  #(*) denotes reverse causation 'm<-y'
  #c1:  '0,0,0' / H1 and H5
  #c2:  '0,0,1' / H2 and H5
  #c3:  '0,1,0' / H1 and H6
  #c4:  '0,1,1' / H2 and H6 - complete mediation
  #c5:  '1,0,0' / H3 and H5 
  #c6:  '1,0,1' / H4 and H5
  #c7:  '1,1,0' / H3 and H6 - colocalization
  #c8:  '1,1,1' / H4 and H6 - partial mediation
  #c9:  '0,0,*' / H1 and H7
  #c10: '0,1,*' / H1 and H8
  #c11: '1,0,*' / H3 and H7
  #c12: '1,1,*' / H3 and H8
  
  preset_odds_index <- return_preset_odds_index()
  output <- posterior_summary(ln_prob_data, ln_prior_c, preset_odds_index)
  colnames(output$ln_post_odds) <- colnames(output$ln_prior_odds) <- colnames(output$ln_post_odds)
  
  #return results
  output$ln_prior_c <- matrix(ln_prior_c, nrow = 1)
  colnames(output$ln_prior_c) <- colnames(output$ln_post_c)
  
  output$ln_prob_data <- ln_prob_data
  output <- output[c("ln_prob_data", "ln_post_c", "ln_post_odds", "ln_prior_c", "ln_prior_odds", "ln_ml")]
  
  if (verbose) { print("Done", quote=F) }
  output
}

## Function to optionally align data (for bmediatR_v0)
align_data_v0 <- function(y, M, X, Z, w,
                          verbose = TRUE) {
  
  # M can have NAs
  overlapping_samples <- Reduce(f = intersect, x = list(names(y), 
                                                        rownames(X),
                                                        rownames(Z), 
                                                        names(w)))
  
  if (length(overlapping_samples) == 0 | !any(overlapping_samples %in% rownames(M))) {
    stop("No samples overlap. Check rownames of M, X, Z and names of y and w.", call. = FALSE)
  } else if (verbose) {
    writeLines(text = c("Number of overlapping samples:", length(overlapping_samples)))
  }
  
  # Return aligned data
  list(y = y[overlapping_samples],
       M = M[overlapping_samples,, drop = FALSE],
       X = X[overlapping_samples,, drop = FALSE],
       Z = Z[overlapping_samples,, drop = FALSE],
       w = w[overlapping_samples])
}

#' Bayesian model selection for mediation analysis function, version 0 
#'
#' This function takes an outcome (y), candidate mediators (M), and a driver as a design matrix (X) to perform a 
#' Bayesian model selection analysis for mediation. Version 0 has more flexibility for specifying the hyper priors
#' on the effect sizes, i.e., the variance explained by the explanatory variables at each step. It also does not
#' have the option for covariates and weights specific to M and y.
#'
#' @param y Vector or single column matrix of an outcome variable. Single outcome variable expected. 
#' Names or rownames must match across M, X, Z, and w (if provided) when align_data = TRUE. If align_data = FALSE,
#' dimensions and order must match across inputs.
#' @param M Vector or matrix of mediator variables. Multiple mediator variables are supported. 
#' Names or rownames must match across y, X, Z and w (if provided) when align_data = TRUE. If align_data = FALSE,
#' dimensions and order must match across inputs.
#' @param X Design matrix of the driver. Names or rownames must match across y, M, Z, and w (if provided) when align_data = TRUE. 
#' If align_data = FALSE, dimensions and order must match across inputs. One common application is for X to represent genetic 
#' information at a QTL, either as founder strain haplotypes or variant genotypes, though X is generalizable to other types of variables.
#' @param Z DEFAULT: NULL. Design matrix of covariates that influence the outcome and mediator variables. 
#' Names or rownames must match across y, M, X, and w (if provided) when align_data = TRUE. If align_data = FALSE,
#' dimensions and order must match across inputs..
#' @param w DEFAULT: NULL. Vector or single column matrix of weights for individuals in analysis. 
#' Names must match across y, M, X, and Z (if provided) when align_data = TRUE. If align_data = FALSE,
#' dimensions and order must match across inputs. A common use would be for an analysis of strain means, where w 
#' is a vector of the number of individuals per strain. If no w is given, observations are equally weighted as 1s.
#' @param tau_sq_mu DEFAULT: 1000. Variance component for the intercept. The DEFAULT represents a diffuse prior, analagous to 
#' a fixed effect term.
#' @param tau_sq_Z DEFAULT: 1000. Variance component for the covariates encoded in Z. The DEFAULT represents a diffuse prior, analagous 
#' to fixed effect terms.
#' @param ln_prior_c DEFAULT: "complete". The prior log case probabilities. See model_info() for description of likelihoods and their
#' combinations into cases. Simplified pre-set options are available, including "complete", "partial", and "reactive".
#' @param align_data DEFAULT: TRUE. If TRUE, expect vector and matrix inputes to have names and rownames, respectively. The overlapping data
#' will then be aligned, allowing the user to not have to reduce data to overlapping samples and order them.
#' @export
#' @examples bmediatR_v0()
bmediatR_v0 <- function(y, M, X, Z = NULL, w = NULL,
                        kappa = rep(0.001, 8),
                        lambda = rep(0.001, 8),
                        tau_sq_mu = rep(1000, 8),
                        tau_sq_Z = rep(1000, 8),
                        phi_sq_X = c(NA,NA,1,0.5,NA,1,NA,0.5),
                        phi_sq_m = c(NA,1,NA,0.5,NA,NA,NA,NA),
                        phi_sq_y = c(NA,NA,NA,NA,NA,NA,1,0.5),
                        ln_prior_c = "complete",
                        options_X = list(sum_to_zero = TRUE, center = FALSE, scale = FALSE),
                        align_data = TRUE,
                        verbose = T) {
  
  #dimension of y
  n <- length(y)
  
  #presets for ln_prior_c; 
  ln_prior_c <- return_ln_prior_c_from_presets(ln_prior_c = ln_prior_c)
  
  #ensure ln_prior_c sum to 1 on probability scale
  if (is.matrix(ln_prior_c)){
    ln_prior_c <- t(apply(ln_prior_c, 1, function(x){x - matrixStats::logSumExp(x)}))
  } else {
    ln_prior_c <- ln_prior_c - matrixStats::logSumExp(ln_prior_c)
  }
  
  #ensure y is a vector
  if (is.matrix(y)) { y <- y[,1] }
  
  #default values for Z
  if (is.null(Z)) { Z <- matrix(NA, n, 0); rownames(Z) <- names(y) }
  
  #ensure X, M, and Z are matrices
  X <- as.matrix(X)
  M <- as.matrix(M)
  Z <- as.matrix(Z)
  
  #default values for w
  if (is.null(w)) { w <- rep(1, n); names(w) <- names(y) }
  
  #optionally align data
  if (align_data) {
    aligned_data <- align_data_v0(y = y, M = M, X = X, Z = Z, w = w,
                                  verbose = verbose)
    y <- aligned_data$y
    M <- aligned_data$M
    X <- aligned_data$X
    Z <- aligned_data$Z
    w <- aligned_data$w
  }

  #drop observations with missing y or X and update n
  complete_y <- !is.na(y)
  complete_X <- !apply(is.na(X), 1, any)
  
  y <- y[complete_y & complete_X]
  M <- M[complete_y & complete_X,,drop=F]
  X <- X[complete_y & complete_X,,drop=F]
  Z <- Z[complete_y & complete_X,,drop=F]
  w <- w[complete_y & complete_X]
  
  n <- length(y)
  
  #drop columns of Z that are invariant
  Z_drop <- which(apply(Z, 2, function(x) var(x)) == 0)
  if (length(Z_drop) > 0) {
    writeLines(paste("Dropping invariants columns from Z:", colnames(Z)[Z_drop]))
    Z <- Z[,-Z_drop]
  }
  
  #dimension of Z
  p <- ncol(Z)
  
  #scale y, M, and Z
  y <- c(scale(y))
  M <- apply(M, 2, scale)
  if (p > 0) { Z <- apply(Z, 2, scale) }
  
  #optionally use sum-to-zero contrast for X
  #recommended when X is a matrix of factors, with a column for every factor level
  if (options_X$sum_to_zero==T){
    C <- sumtozero_contrast(ncol(X))
    X <- X%*%C
  }
  
  #optionally center and scale X
  X <- apply(X, 2, scale, center = options_X$center, scale = options_X$scale)
  
  #dimension of X
  d <- ncol(X)
  
  #column design matrix for mu
  ones <- matrix(1, n)
  
  #begin Bayesian calculations
  if (verbose) { print("Initializing", quote = FALSE) }
  
  #identify likelihoods that are not supported by the prior
  #will not compute cholesky or likelihood for these
  calc_ln_prob_data <- rep(NA, 8) 
  calc_ln_prob_data[1] <- any(!is.infinite(ln_prior_c[c(1,3,9,10)]))
  calc_ln_prob_data[2] <- any(!is.infinite(ln_prior_c[c(2,4)]))
  calc_ln_prob_data[3] <- any(!is.infinite(ln_prior_c[c(5,7,11,12)]))
  calc_ln_prob_data[4] <- any(!is.infinite(ln_prior_c[c(6,8)]))
  calc_ln_prob_data[5] <- any(!is.infinite(ln_prior_c[c(1,2,5,6)]))
  calc_ln_prob_data[6] <- any(!is.infinite(ln_prior_c[c(3,4,7,8)]))
  calc_ln_prob_data[7] <- any(!is.infinite(ln_prior_c[c(9,11)]))
  calc_ln_prob_data[8] <- any(!is.infinite(ln_prior_c[c(10,12)]))
  
  #likelihood models for all hypothesis
  #hypotheses encoded by presence (1) or absence (0) of 'X->y, X->m, m->y' edges on the DAG
  #(*) denotes reverse causation 'm<-y'
  #H1: '0,-,0' / y does not depend on X or m
  #H2: '0,-,1' / y depends on m but not X
  #H3: '1,-,0' / y depends on X but not m
  #H4: '1,-,1' / y depends on X and m
  #H5: '-,0,-' / m does not depend on X
  #H6: '-,1,-' / m depends on X
  #H7: '-,0,*' / m depends on y but not X
  #H8: '-,1,*' / m depends on X and y
  #all include covariates Z
  
  #design matrices for H1,H3,H5-H8 complete cases (do not depend on m)
  X1 <- cbind(ones, Z)
  X3 <- cbind(ones, X, Z)
  X5 <- X1
  X6 <- X3
  X7 <- cbind(X1, y)
  X8 <- cbind(X3, y)
  
  #check if all scale hyperparameters are identical for H1 and H5
  #implies sigma1 and sigma5 identical, used to reduce computations
  sigma5_equal_sigma1 <- all(lambda[1]==lambda[5],
                             tau_sq_mu[1] == tau_sq_mu[5], 
                             tau_sq_Z[1] == tau_sq_Z[5])
  
  #check if all scale hyperparameters are identical for H3 and H6
  #implies sigma3 and sigma6 identical, used to reduce computations
  sigma6_equal_sigma3 <- all(lambda[3]==lambda[6],
                             tau_sq_mu[3] == tau_sq_mu[6], 
                             tau_sq_Z[3] == tau_sq_Z[6])
  
  #prior variance matrices (diagonal) for H1-H8 
  v1 <- c(tau_sq_mu[1], rep(tau_sq_Z[1], p))
  v2 <- c(tau_sq_mu[2], rep(tau_sq_Z[2], p), phi_sq_m[2])
  v3 <- c(tau_sq_mu[3], rep(phi_sq_X[3], d), rep(tau_sq_Z[3], p))
  v4 <- c(tau_sq_mu[4], rep(phi_sq_X[4], d), rep(tau_sq_Z[4], p), phi_sq_m[4])
  v7 <- c(tau_sq_mu[7], rep(tau_sq_Z[7], p), phi_sq_y[7])
  v8 <- c(tau_sq_mu[8], rep(phi_sq_X[8], d), rep(tau_sq_Z[8], p), phi_sq_y[8])
  
  if (!sigma5_equal_sigma1 | !calc_ln_prob_data[1]){
    v5 <- c(tau_sq_mu[5], rep(tau_sq_Z[5], p))
  }
  
  if (!sigma6_equal_sigma3 | !calc_ln_prob_data[3]){
    v6 <- c(tau_sq_mu[6], rep(phi_sq_X[6], d), rep(tau_sq_Z[6], p))
  }
  
  #scale matrices for H1,H3,H5-H8 complete cases (do not depend on m)
  sigma1 <- crossprod(sqrt(lambda[1]*v1)*t(X1))
  sigma3 <- crossprod(sqrt(lambda[3]*v3)*t(X3))
  sigma7 <- crossprod(sqrt(lambda[7]*v7)*t(X7))
  sigma8 <- crossprod(sqrt(lambda[8]*v8)*t(X8))
  
  diag(sigma1) <- diag(sigma1) + lambda[1]/w
  diag(sigma3) <- diag(sigma3) + lambda[3]/w
  diag(sigma7) <- diag(sigma7) + lambda[7]/w
  diag(sigma8) <- diag(sigma8) + lambda[8]/w
  
  if (!sigma5_equal_sigma1 | !calc_ln_prob_data[1]){
    sigma5 <- crossprod(sqrt(lambda[5]*v5)*t(X5))
    diag(sigma5) <- diag(sigma5) + lambda[5]/w
  }
  
  if (!sigma6_equal_sigma3 | !calc_ln_prob_data[3]){
    sigma6 <- crossprod(sqrt(lambda[6]*v6)*t(X6))
    diag(sigma6) <- diag(sigma6) + lambda[6]/w
  }
  
  #object to store likelihoods
  ln_prob_data=matrix(-Inf, ncol(M), 8)
  rownames(ln_prob_data) <- colnames(M)
  colnames(ln_prob_data) <- c("0,-,0", 
                              "0,-,1", 
                              "1,-,0", 
                              "1,-,1",
                              "-,0,-", 
                              "-,1,-",
                              '-,0,*',
                              '-,1,*')
  
  #identify batches of M that have the same pattern of missing values
  missing_m <- batch_cols(M)
  
  #iterate over batches of M with same pattern of missing values
  if(verbose){print("Iterating", quote=F)}
  counter <- 0
  
  for (b in 1:length(missing_m)){
    #subset to non-missing observations
    index <- rep(T, length(y))
    index[missing_m[[b]]$omit] <- F
    
    if (any(index)){
      y_subset <- y[index]
      w_subset <- w[index]
      
      #cholesky matrices for H1,H3,H5-H8 non-missing observations (do not depend on m)
      if (calc_ln_prob_data[1]){sigma1_chol_subset <- chol(sigma1[index,index])}
      if (calc_ln_prob_data[3]){sigma3_chol_subset <- chol(sigma3[index,index])}
      if (calc_ln_prob_data[7]){sigma7_chol_subset <- chol(sigma7[index,index])}
      if (calc_ln_prob_data[8]){sigma8_chol_subset <- chol(sigma8[index,index])}
      
      if (sigma5_equal_sigma1 & calc_ln_prob_data[1]){
        sigma5_chol_subset <- sigma1_chol_subset
      } else if (calc_ln_prob_data[5]){
        sigma5_chol_subset <- chol(sigma5[index,index])
      }
      
      if (sigma6_equal_sigma3 & calc_ln_prob_data[3]){
        sigma6_chol_subset <- sigma3_chol_subset
      } else if (calc_ln_prob_data[6]){
        sigma6_chol_subset <- chol(sigma6[index,index])
      }
      
      #compute H1 and H3 outside of the mediator loop (invariant)
      if (calc_ln_prob_data[1]){ln_prob_data1 <- dmvt_chol(y_subset, sigma_chol=sigma1_chol_subset, df = kappa[1])}
      if (calc_ln_prob_data[3]){ln_prob_data3 <- dmvt_chol(y_subset, sigma_chol=sigma3_chol_subset, df = kappa[3])}
      
      #iterate over mediators
      for (i in missing_m[[b]]$cols){
        counter <- counter + 1
        if (counter%%1000==0 & verbose){print(paste(counter, "of", ncol(M)), quote=F)}
        
        #set current mediator non-missing observations
        m_subset <- M[index,i]
        
        #design matrix for H2 and H4 non-missing observations
        X2_subset <- cbind(X1[index,,drop=F], m_subset)
        X4_subset <- cbind(X3[index,,drop=F], m_subset)
        
        #scale and cholesky matrices for H2 and H4 non-missing observations
        sigma2_subset <- crossprod(sqrt(lambda[2]*v2)*t(X2_subset))
        sigma4_subset <- crossprod(sqrt(lambda[4]*v4)*t(X4_subset))
        
        diag(sigma2_subset) <- lambda[2]/w_subset + diag(sigma2_subset)
        diag(sigma4_subset) <- lambda[4]/w_subset + diag(sigma4_subset)
        
        if (calc_ln_prob_data[2]){sigma2_chol_subset <- chol(sigma2_subset)}
        if (calc_ln_prob_data[4]){sigma4_chol_subset <- chol(sigma4_subset)}
        
        #compute likelihoods for H1-H8
        if (calc_ln_prob_data[1]){ln_prob_data[i,1] <- ln_prob_data1}
        if (calc_ln_prob_data[3]){ln_prob_data[i,3] <- ln_prob_data3}
        
        if (calc_ln_prob_data[2]){ln_prob_data[i,2] <- dmvt_chol(y_subset, sigma_chol=sigma2_chol_subset, df = kappa[2])}
        if (calc_ln_prob_data[4]){ln_prob_data[i,4] <- dmvt_chol(y_subset, sigma_chol=sigma4_chol_subset, df = kappa[4])}
        if (calc_ln_prob_data[5]){ln_prob_data[i,5] <- dmvt_chol(m_subset, sigma_chol=sigma5_chol_subset, df = kappa[5])}
        if (calc_ln_prob_data[6]){ln_prob_data[i,6] <- dmvt_chol(m_subset, sigma_chol=sigma6_chol_subset, df = kappa[6])}
        if (calc_ln_prob_data[7]){ln_prob_data[i,7] <- dmvt_chol(m_subset, sigma_chol=sigma7_chol_subset, df = kappa[7])}
        if (calc_ln_prob_data[8]){ln_prob_data[i,8] <- dmvt_chol(m_subset, sigma_chol=sigma8_chol_subset, df = kappa[8])}
      }
    }
  }
  
  #compute posterior probabilities for all cases
  #compute posterior odds for specified combinations of cases
  #cases encoded by presence (1) or absence (0) of 'X->y, X->m, m->y' edges on the DAG
  #(*) denotes reverse causation 'm<-y'
  #c1:  '0,0,0' / H1 and H5
  #c2:  '0,0,1' / H2 and H5
  #c3:  '0,1,0' / H1 and H6
  #c4:  '0,1,1' / H2 and H6 - complete mediation
  #c5:  '1,0,0' / H3 and H5 
  #c6:  '1,0,1' / H4 and H5
  #c7:  '1,1,0' / H3 and H6 - colocalization
  #c8:  '1,1,1' / H4 and H6 - partial mediation
  #c9:  '0,0,*' / H1 and H7
  #c10: '0,1,*' / H1 and H8
  #c11: '1,0,*' / H3 and H7
  #c12: '1,1,*' / H3 and H8
  
  preset_odds_index <- return_preset_odds_index()
  output <- posterior_summary(ln_prob_data, ln_prior_c, preset_odds_index)
  colnames(output$ln_post_odds) <- colnames(output$ln_prior_odds) <- colnames(output$ln_post_odds)
  
  #return results
  output$ln_prior_c <- matrix(ln_prior_c, nrow = 1)
  colnames(output$ln_prior_c) <- colnames(output$ln_post_c)
    
  output$ln_prob_data <- ln_prob_data
  output <- output[c("ln_prob_data", "ln_post_c", "ln_post_odds", "ln_prior_c", "ln_prior_odds", "ln_ml")]
  
  if (verbose) {print("Done", quote=F)}
  output
}

#' Definions and encodings for models included in the mediation analysis through Bayesian model selection 
#'
#' This function prints out a table describing the models considered in the mediation analysis.
#'
#' @export
#' @examples model_info()
model_info <- function(){
  writeLines(c("likelihood models for all hypotheses",
               "hypotheses encoded by presence (1) or absence (0) of 'X->y, X->m, m->y' edges on the DAG",
               "(*) denotes reverse causation 'm<-y'",
               "H1: '0,-,0' / y does not depend on X or m",
               "H2: '0,-,1' / y depends on m but not X",
               "H3: '1,-,0' / y depends on X but not m",
               "H4: '1,-,1' / y depends on X and m",
               "H5: '-,0,-' / m does not depend on X",
               "H6: '-,1,-' / m depends on X",
               "H7: '-,0,*' / m depends on y but not X",
               "H8: '-,1,*' / m depends on X and y",
               "",
               "combinations of hypotheses for all cases",
               "cases encoded by presence (1) or absence (0) of 'X->y, X->m, m->y' edges on the DAG",
               "(*) denotes reverse causation 'm<-y'",
               "c1:  '0,0,0' / H1 and H5",
               "c2:  '0,0,1' / H2 and H5",
               "c3:  '0,1,0' / H1 and H6",
               "c4:  '0,1,1' / H2 and H6 - complete mediation",
               "c5:  '1,0,0' / H3 and H5 ",
               "c6:  '1,0,1' / H4 and H5",
               "c7:  '1,1,0' / H3 and H6 - colocalization",
               "c8:  '1,1,1' / H4 and H6 - partial mediation",
               "c9:  '0,0,*' / H1 and H7",
               "c10: '0,1,*' / H1 and H8",
               "c11: '1,0,*' / H3 and H7",
               "c12: '1,1,*' / H3 and H8"))
}

#' Calculate empirical priors based on the relationships in the data between X, M, and y
#'
#' This function estimates log prior model probabilities that maximize the marginal likelihoods
#' of the data from a prior run of bmediatR(). 
#' 
#' @param ln_prob_data Marginal likelihoods from bmediatR fit.
#' @param model DEFAULT: "complete". If "complete", the marginal likelihoods evaluating the 
#' effects of X on M and Y. If "partial" is specified, the effects of X on M and Y are fixed
#' as present. A "reactive" mode is currently not supported because of the challenge in
#' distinguishing the directionality of a relationship.
#' @export
#' @examples model_info()
estimate_empirical_prior <- function(ln_prob_data,
                                     model = c("complete", "partial")) {
  model <- model[1]
  
  if (model == "complete") {
    marginal_likelihood <- function(x){
      x <- -log(1+exp(-x))
      
      ln_prior_c <- rep(-Inf, 12)
      ln_prior_c[1] <- VGAM::log1mexp(-x[1]) + VGAM::log1mexp(-x[2]) + VGAM::log1mexp(-x[3])
      ln_prior_c[2] <- VGAM::log1mexp(-x[1]) + VGAM::log1mexp(-x[2]) + x[3]
      ln_prior_c[3] <- VGAM::log1mexp(-x[1]) + x[2] + VGAM::log1mexp(-x[3])
      ln_prior_c[4] <- VGAM::log1mexp(-x[1]) + x[2] + x[3]
      ln_prior_c[5] <- x[1] + VGAM::log1mexp(-x[2]) + VGAM::log1mexp(-x[3])
      ln_prior_c[6] <- x[1] + VGAM::log1mexp(-x[2]) + x[3]
      ln_prior_c[7] <- x[1] + x[2] + VGAM::log1mexp(-x[3])
      ln_prior_c[8] <- x[1] + x[2] + x[3]
      
      -sum(posterior_summary(ln_prob_data, ln_prior_c, list(1))$ln_ml)
    }
    
    empirical_prior <- optim(c(0,0,0), marginal_likelihood)
    
    x <- -log(1+exp(-empirical_prior$par))
    
    ln_prior_c <- rep(-Inf, 12)
    ln_prior_c[1] <- VGAM::log1mexp(-x[1]) + VGAM::log1mexp(-x[2]) + VGAM::log1mexp(-x[3])
    ln_prior_c[2] <- VGAM::log1mexp(-x[1]) + VGAM::log1mexp(-x[2]) + x[3]
    ln_prior_c[3] <- VGAM::log1mexp(-x[1]) + x[2] + VGAM::log1mexp(-x[3])
    ln_prior_c[4] <- VGAM::log1mexp(-x[1]) + x[2] + x[3]
    ln_prior_c[5] <- x[1] + VGAM::log1mexp(-x[2]) + VGAM::log1mexp(-x[3])
    ln_prior_c[6] <- x[1] + VGAM::log1mexp(-x[2]) + x[3]
    ln_prior_c[7] <- x[1] + x[2] + VGAM::log1mexp(-x[3])
    ln_prior_c[8] <- x[1] + x[2] + x[3]
  }
  else if (model == "partial") {
    marginal_likelihood <- function(x){
      x <- -log(1+exp(-x))
      
      ln_prior_c <- rep(-Inf, 12)
      ln_prior_c[5] <- VGAM::log1mexp(-x[1]) + VGAM::log1mexp(-x[2])
      ln_prior_c[6] <- VGAM::log1mexp(-x[1]) + x[2]
      ln_prior_c[7] <- x[1] + VGAM::log1mexp(-x[2])
      ln_prior_c[8] <- x[1] + x[2]
      
      -sum(posterior_summary(ln_prob_data, ln_prior_c, list(1))$ln_ml)
    }
    
    empirical_prior <- optim(c(0,0), marginal_likelihood)
    
    x <- -log(1+exp(-empirical_prior$par))
    
    ln_prior_c <- rep(-Inf, 12)
    ln_prior_c[5] <- VGAM::log1mexp(-x[1]) + VGAM::log1mexp(-x[2])
    ln_prior_c[6] <- VGAM::log1mexp(-x[1]) + x[2]
    ln_prior_c[7] <- x[1] + VGAM::log1mexp(-x[2])
    ln_prior_c[8] <- x[1] + x[2]
  }
  
  ln_prior_c <- matrix(ln_prior_c, nrow = 1)
  colnames(ln_prior_c) <- c("0,0,0", "0,0,1", "0,1,0", "0,1,1", 
                            "1,0,0", "1,0,1", "1,1,0", "1,1,1", 
                            "0,0,*", "0,1,*", "1,0,*", "1,1,*")
  ln_prior_c
}

