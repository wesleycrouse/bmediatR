get_approx_pval <- function(med_bf_object,
                            annot) {
  
  bf_dat <- data.frame(protein.id = names(med_bf_object$lnBF), lnBF = med_bf_object$lnBF) %>%
    left_join(annot)
  
  approx_pval <- nrow(bf_dat)
  for (i in 1:nrow(bf_dat)) {
    approx_pval[i] <- mean(c(bf_dat$lnBF[i] < (bf_dat %>% filter(chr != bf_dat$chr[i]) %>% pull(lnBF)), TRUE))
  }
  bf_dat$approx_pval <- approx_pval
  bf_dat
}

get_perm_pval <- function(y, 
                          M, 
                          X, 
                          Z = NULL,
                          num_perm = 1000,
                          verbose = FALSE,
                          ...) {
  
  actual_bf <- mediation_bf(y = y, 
                            M = M, 
                            X = X,
                            Z = Z,
                            verbose = verbose,
                            ...)
  
  if (is.null(Z)) { Z <- matrix(1, nrow = nrow(M)); rownames(Z) <- rownames(M) }
  
  within_cov <- Z %>% unique
  original_name <- NULL
  perm_mat <- NULL
  for (i in 1:nrow(within_cov)) {
    sub_Z <- Z[apply(Z, 1, function(x) paste(x, collapse = "")) == paste(within_cov[i,], collapse = ""),]
    sub_original_name <- rownames(sub_Z)
    sub_perm_mat <- matrix(NA, nrow = nrow(sub_Z), ncol = num_perm)
    for (j in 1:num_perm) {
      sub_perm_mat[,j] <- sample(sub_original_name)
    }
    original_name <- c(original_name, sub_original_name)
    perm_mat <- rbind(perm_mat, sub_perm_mat)
  }
  
  rownames(perm_mat) <- original_name
  perm_mat <- perm_mat[rownames(M),]
  
  perm_pval <- rep(NA, ncol(M))
  perm_bf_mat <- matrix(NA, nrow = num_perm, ncol = ncol(M))
  colnames(perm_bf_mat) <- colnames(M)
  for(i in 1:num_perm) {
    ## Rename things
    perm_M <- M
    rownames(perm_M) <- as.character(perm_mat[,i])
    for (j in 1:ncol(M)) {
      perm_med <- mediation_bf(y = y, 
                               M = perm_M[,j,drop = FALSE][rownames(M),], 
                               X = X,
                               Z = Z, 
                               verbose = verbose)
      perm_bf_mat[i, j] <- perm_med$lnBF
    } 
    print(paste("Perm", i, "done out of", num_perm))
  }
  
  perm_pval <- sapply(1:ncol(perm_bf_mat), function(i) mean(actual_bf$lnBF[i] < perm_bf_mat[,i]))
  bf_dat <- data.frame(protein.id = names(actual_bf$lnBF), lnBF = actual_bf$lnBF, perm_pval = perm_pval)
  bf_dat
}
