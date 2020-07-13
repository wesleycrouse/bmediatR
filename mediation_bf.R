require(matrixStats)

sumtozero_contrast <- function(k){
  u <- 1/((k-1)^(0.5))
  v <- (-1+(k^(0.5)))*((k-1)^(-1.5))
  w <- (k-2)*v+u
  C <- matrix(-v, k-1, k-1)
  diag(C) <- rep(w, k-1)
  rbind(C, rep(-u, k-1))
}

#adapted from mvtnorm::dmvt
dmvt_chol <- function (x, sigma_chol, df){
  if (is.vector(x)){x <- matrix(x, nrow=length(x))}
  n <- nrow(x)
  dec <- sigma_chol
  R_x_m <- backsolve(dec, x, transpose = TRUE)
  rss <- colSums(R_x_m^2)
  lgamma(0.5*(df+n)) - (lgamma(0.5*df) + sum(log(diag(dec))) + 0.5*n*log(pi*df)) - 0.5*(df+n)*log1p(rss/df)
}

#adapted from qtl2::batch_cols
batch_cols <- function(mat){
  mat <- !is.finite(mat)
  n <- nrow(mat)
  all_true <- rep(TRUE, n)
  result <- NULL
  n_na <- colSums(mat)
  no_na <- (n_na == 0)
  if (any(no_na)){result <- list(list(cols=which(no_na), omit=numeric(0)))}
  one_na <- (n_na == 1)
  if (any(one_na)){
    wh <- apply(mat[,one_na, drop = FALSE], 2, which)
    spl <- split(which(one_na), wh)
    part2 <- lapply(seq_along(spl), function(i){list(cols=as.numeric(spl[[i]]), omit=as.numeric(names(spl)[i]))})
    if (is.null(result)){
      result <- part2
    } else {
      result <- c(result, part2)
    }
  }
  other_cols <- !(no_na | one_na)
  if (any(other_cols)){
    other_cols <- seq_len(ncol(mat))[other_cols]
    pat <- apply(mat[,other_cols, drop = FALSE], 2, function(a) paste(which(a), collapse = ":"))
    u <- unique(pat)
    part3 <- lapply(u, function(a) list(cols=other_cols[pat==a], omit=as.numeric(strsplit(a, ":")[[1]])))
  } else {
    part3 <- NULL
  }
  if (is.null(result)){
    part3
  } else {
    c(result, part3)
  }
}

mediation_bf <- function(y, M, X, Z = NULL, w = NULL,
                         kappa = rep(0.001, 4),
                         lambda = rep(0.001, 4),
                         tau_sq_mu = rep(1000, 4),
                         tau_sq_Z = rep(1000, 4),
                         phi_sq_X = c(NA, 1, 1, 0.5),
                         phi_sq_m = 0.5,
                         p_c = c(1/3, 1/3, 1/3),
                         p_med = 0.5,
                         verbose = T){
  if(verbose){print("Initializing", quote=F)}
  
  #dimensions
  n <- nrow(X)
  d <- ncol(X)-1
  p <- ncol(Z)
  
  #default values for w and Z
  if (is.null(w)){w <- rep(1, n)}
  if (is.null(Z)){Z <- matrix(NA, n, 0); p <- 0}
  
  #ensure M and Z are matrices
  M <- as.matrix(M)
  Z <- as.matrix(Z)
  
  #drop observations with missing y
  complete.y <- !is.na(y)
  
  y <- y[complete.y]
  M <- M[complete.y,,drop=F]
  X <- X[complete.y,,drop=F]
  Z <- Z[complete.y,,drop=F]
  w <- w[complete.y]
  
  #scale y, M, and Z
  y <- c(scale(y))
  M <- apply(M, 2, scale)
  if (p > 0){Z <- apply(Z, 2, scale)}
  
  #precomputed quantities for design matrices
  C <- sumtozero_contrast(ncol(X))
  XC <- X%*%C
  ones <- matrix(1, nrow(XC))
  
  #design matrices for H1-H3 complete case
  #H1: mediator m does not depend on X
  #H2: trait y depends on X but not on mediator m
  #H3: mediator m depends on X
  #H4: trait y depends on X and mediator m
  #all have covariates Z
  X1 <- cbind(ones, Z)
  X2 <- cbind(ones, XC, Z)
  X3 <- X2
  
  #check if all scale hyperparameters are identical for H2 and H3
  #implies sigma2 and sigma3 identical, used to reduce computations
  sigma3_equal_sigma2 <- all(lambda[2]==lambda[3],
                             tau_sq_mu[2] == tau_sq_mu[3], 
                             phi_sq_X[2] == phi_sq_X[3],
                             tau_sq_Z[2] == tau_sq_Z[3])
  
  #prior variance matrices (diagonal) for H1-H4 
  v1 <- c(tau_sq_mu[1], rep(tau_sq_Z[1], p))
  v2 <- c(tau_sq_mu[2], rep(phi_sq_X[2], d), rep(tau_sq_Z[2], p))
  v4 <- c(tau_sq_mu[4], rep(phi_sq_X[4], d), rep(tau_sq_Z[4], p), phi_sq_m)
  
  if (!sigma3_equal_sigma2){
    v3 <- c(tau_sq_mu[3], rep(phi_sq_X[3], d), rep(tau_sq_Z[3], p))
  }
  
  #scale matrices for H1-H3 complete case
  sigma1 <- crossprod(sqrt(lambda[1]*v1)*t(X1))
  sigma2 <- crossprod(sqrt(lambda[2]*v2)*t(X2))
  
  diag(sigma1) <- diag(sigma1) + lambda[1]/w
  diag(sigma2) <- diag(sigma2) + lambda[2]/w
  
  if (!sigma3_equal_sigma2){
    sigma3 <- crossprod(sqrt(lambda[3]*v3)*t(X3))
    diag(sigma3) <- diag(sigma3) + lambda[3]/w
  }
  
  #object to store likelihoods
  lnp_data_H=matrix(NA, ncol(M), 4)
  rownames(lnp_data_H) <- colnames(M)
  colnames(lnp_data_H) <- c("H1", "H2", "H3", "H4")
  
  #identify batches of M that have the same pattern of missing values
  missing_m <- batch_cols(M)
  
  #iterate over batches of M with same pattern of missing values
  if(verbose){print("Iterating", quote=F)}
  counter <- 0
  
  for (b in 1:length(missing_m)){
    #subset to non-missing observations
    index <- rep(T, length(y))
    index[missing_m[[b]]$omit] <- F
    
    y_subset <- y[index]
    w_subset <- w[index]
    
    #cholesky matrices for H1-H3 non-missing observations
    sigma1_chol_subset <- chol(sigma1[index,index])
    sigma2_chol_subset <- chol(sigma2[index,index])
    
    if (sigma3_equal_sigma2){
      sigma3_chol_subset <- sigma2_chol_subset
    } else {
      sigma3_chol_subset <- chol(sigma3[index,index])
    }
    
    #compute H2 outside of the mediator loop (invariant)
    lnp_data_H2 <- dmvt_chol(y_subset, sigma_chol=sigma2_chol_subset, df = kappa[2])
    
    #iterate over mediators
    for (i in missing_m[[b]]$cols){
      counter <- counter + 1
      if (counter%%1000==0 & verbose){print(paste(counter, "of", ncol(M)), quote=F)}
      
      #set current mediator non-missing observations
      m_subset <- M[index,i]
      
      #design matrix for H4 non-missing observations
      X4_subset <- cbind(X2[index,,drop=F], m_subset)
      
      #scale and cholesky matrices for H4 non-missing observations
      sigma4_subset <- crossprod(sqrt(lambda[4]*v4)*t(X4_subset))
      diag(sigma4_subset) <- lambda[4]/w_subset + diag(sigma4_subset)
      sigma4_chol_subset <- chol(sigma4_subset)
      
      #compute likelihoods for H1-H4
      lnp_data_H[i,2] <- lnp_data_H2
      lnp_data_H[i,1] <- dmvt_chol(m_subset, sigma_chol=sigma1_chol_subset, df = kappa[1])
      lnp_data_H[i,3] <- dmvt_chol(m_subset, sigma_chol=sigma3_chol_subset, df = kappa[3])
      lnp_data_H[i,4] <- dmvt_chol(y_subset, sigma_chol=sigma4_chol_subset, df = kappa[4])
    }
  }
  
  #compute Bayes factors
  lnBF_numerator_med <- lnp_data_H[,3] + lnp_data_H[,4]
  
  lnBF_denominator_med <- rbind(lnp_data_H[,1] + lnp_data_H[,2] + log(p_c[1]),
                                lnp_data_H[,1] + lnp_data_H[,4] + log(p_c[2]),
                                lnp_data_H[,3] + lnp_data_H[,2] + log(p_c[3]))
  lnBF_denominator_med <- apply(lnBF_denominator_med, 2, matrixStats::logSumExp)
  
  lnBF_med<- lnBF_numerator_med - lnBF_denominator_med
  
  #computer posterior probabilities of mediation
  ln_post_med <- lnBF_med + log(p_med) - apply(rbind(lnBF_med+log(p_med), log(1-p_med)) , 2, matrixStats::logSumExp)
  
  #compute co-local Bayes factors
  lnBF_numerator_coloc <- lnp_data_H[,3] + lnp_data_H[,2]
  
  lnBF_denominator_coloc <- rbind(lnp_data_H[,1] + lnp_data_H[,2] + log(p_c[1]),
                                  lnp_data_H[,1] + lnp_data_H[,4] + log(p_c[2]),
                                  lnp_data_H[,3] + lnp_data_H[,4] + log(p_c[3]))
  lnBF_denominator_coloc <- apply(lnBF_denominator_coloc, 2, matrixStats::logSumExp)
  
  lnBF_coloc <- lnBF_numerator_coloc - lnBF_denominator_coloc
  
  if(verbose){print("Done", quote=F)}
  list(lnBF_med=lnBF_med, lnBF_coloc=lnBF_coloc, lnp_data_H=lnp_data_H, ln_post_med=ln_post_med)
}







