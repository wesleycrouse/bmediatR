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

#' Summaries of posterior model probabilities function 
#'
#' This function takes the log posterior probability of the data (posterior likelihood) for the various models, the log prior model probabilities, and
#' returns log posterior odds
#'
#' @param ln_prob_data Log posterior likelihoods under the various models, returned by bmediatR().
#' @param ln_prior_c Log prior model probabilities. If posterior_summary() is being used for a non-default posterior odds
#' summary, the log prior model probabilities used with bmediatR() are stored in its output.
#' @param c_numerator The index of models to be summed in the numberator of the posterior odds. Models and their order provided with
#' model_summary().
#' @export
#' @examples bmediatR()
posterior_summary <- function(ln_prob_data, ln_prior_c, c_numerator){
  #function to compute log odds from log probabilities
  ln_odds <- function(ln_p, numerator){
    ln_odds_numerator <- apply(ln_p[,numerator,drop=F], 1, matrixStats::logSumExp)
    ln_odds_denominator <- apply(ln_p[,-numerator,drop=F], 1, matrixStats::logSumExp)
    ln_odds <- ln_odds_numerator -ln_odds_denominator
  }
  
  #ensure c_numerator is a list
  if (!is.list(c_numerator)){
    c_numerator <- list(c_numerator)
  }
  
  #presets for ln_prior_c; 
  if (ln_prior_c[1]=="complete"){
    ln_prior_c <- c(rep(0,8), rep(-Inf,4))
  } else if (ln_prior_c[1]=="partial"){
    ln_prior_c <- c(rep(-Inf,4), rep(0,4), rep(-Inf,4))
  } else if (ln_prior_c[1]=="reactive"){
    ln_prior_c <- rep(0,12)
  }
  
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
  ln_post_c <- ln_post_c - apply(ln_post_c, 1, matrixStats::logSumExp)
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
  colnames(ln_prior_odds) <- c_numerator
  
  #compute posterior odds for each combination of cases
  ln_post_odds <- sapply(c_numerator, ln_odds, ln_p=ln_post_c)
  ln_post_odds <- matrix(ln_post_odds, ncol=length(c_numerator))
  rownames(ln_post_odds) <- rownames(ln_post_c)
  colnames(ln_post_odds) <- c_numerator
  
  #return results
  list(ln_post_c=ln_post_c, ln_post_odds=ln_post_odds, ln_prior_odds=ln_prior_odds)
}

#' Bayesian model selection for mediation analysis function 
#'
#' This function takes an outcome, mediator(s), and a driver as a design matrix to perform a Bayesian model selection analysis for mediation.
#'
#' @param y Vector or single column matrix of an outcome variable. Single outcome variable expected. 
#' Names or rownames must match across M and X, and Z and w (if provided).
#' @param M Vector or matrix of mediator variables. Multiple mediator variables are supported. 
#' Names or rownames must match across y and X, and Z and w (if provided).
#' @param X Design matrix of the driver. Names or rownames must match across y and M, and Z and w (if provided).
#' One common application is for X to represent genetic information at a QTL as either founder strain haplotypes
#' in a multiparental population or variant genotypes, though it is generalizable to types of variables.
#' @param Z DEFAULT: NULL. Design matrix of covariates that influence the outcome variable. Names or rownames must match across y, M and X, 
#' and w (if provided).
#' @param w DEFAULT: NULL. Vector of weights for individuals in analysis. Names or rownames must match across y, M X, and Z (if provided).
#' A common use would be for an analysis of strain means, where w would be a vector of the number of individuals per strain.
#' If no w is given, observations are equally weighted as 1.
#' @export
#' @examples bmediatR()
bmediatR <- function(y, M, X, Z = NULL, w = NULL,
                     kappa = rep(0.001, 8),
                     lambda = rep(0.001, 8),
                     tau_sq_mu = rep(1000, 8),
                     tau_sq_Z = rep(1000, 8),
                     phi_sq_X = c(NA,NA,1,0.5,NA,1,NA,0.5),
                     phi_sq_m = c(NA,1,NA,0.5,NA,NA,NA,NA),
                     phi_sq_y = c(NA,NA,NA,NA,NA,NA,1,0.5),
                     ln_prior_c = "complete",
                     options_X = list(sum_to_zero=T, center=F, scale=F),
                     verbose = T) {
  if(verbose){print("Initializing", quote=F)}
  
  #dimensions
  n <- length(y)
  d <- ncol(X)
  p <- ncol(Z)
  
  #presets for ln_prior_c; 
  if (ln_prior_c=="complete"){
    ln_prior_c <- c(rep(0,8), rep(-Inf,4))
  } else if (ln_prior_c=="partial"){
    ln_prior_c <- c(rep(-Inf,4), rep(0,4), rep(-Inf,4))
  } else if (ln_prior_c=="reactive"){
    ln_prior_c <- rep(0,12)
  }
  
  #ensure ln_prior_c sum to 1 on probability scale
  if (is.matrix(ln_prior_c)){
    ln_prior_c <- t(apply(ln_prior_c, 1, function(x){x - matrixStats::logSumExp(x)}))
  } else {
    ln_prior_c <- ln_prior_c - matrixStats::logSumExp(ln_prior_c)
  }
  
  #default values for w and Z
  if (is.null(w)){w <- rep(1, n)}
  if (is.null(Z)){Z <- matrix(NA, n, 0); p <- 0}
  
  #ensure X, M, and Z are matrices
  X <- as.matrix(X)
  M <- as.matrix(M)
  Z <- as.matrix(Z)
  
  #drop observations with missing y or X and update n
  complete_y <- !is.na(y)
  complete_X <- !apply(is.na(X), 1, any)
  
  y <- y[complete_y & complete_X]
  M <- M[complete_y & complete_X,,drop=F]
  X <- X[complete_y & complete_X,,drop=F]
  Z <- Z[complete_y & complete_X,,drop=F]
  w <- w[complete_y & complete_X]
  
  n <- length(y)
  
  #scale y, M, and Z
  y <- c(scale(y))
  M <- apply(M, 2, scale)
  if (p > 0){Z <- apply(Z, 2, scale)}
  
  #optionally use sum-to-zero contrast for X
  #recommended when X is a matrix of factors, with a column for every factor level
  if (options_X$sum_to_zero==T){
    C <- sumtozero_contrast(ncol(X))
    X <- X%*%C
    d <- ncol(X)
  }
  
  #optionally center and scale X
  X <- apply(X, 2, scale, center=options_X$center, scale=options_X$scale)
  
  #column design matrix for mu
  ones <- matrix(1, n)
  
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
  
  output <- posterior_summary(ln_prob_data, ln_prior_c, list(c(4,8), 8, 4, 7, 9:12, c(4:8,11,12)))
  colnames(output$ln_post_odds) <- c("mediation", "partial", "complete", "colocal", "reactive", "y_depends_x")
  colnames(output$ln_prior_odds) <- colnames(output$ln_post_odds)
  
  #return results
  output$ln_prior_c <- matrix(ln_prior_c, nrow = 1)
  colnames(output$ln_prior_c) <- colnames(output$ln_post_c)
    
  output$ln_prob_data <- ln_prob_data
  output <- output[c(5,1,2,4,3)]
  
  if (verbose) {print("Done", quote=F)}
  output
}

#' Definions and encodings for models included in the mediation analysis through Bayesian model selection 
#'
#' This function prints out a table describing the models considered in the mediation analysis
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