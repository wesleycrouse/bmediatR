#' Function that calculate approximate p-values for Bayes factors
#'
#' This function takes Bayes factor results from mediation_bf() and calculates approximate p-values based
#' on comparisons to genome results, using the assumption that the vast majority of candidates are not
#' mediators and can thus be used to characterize a null distribution of Bayes factors.
#'
#' @export
#' @examples get_approx_pval()
get_approx_pval <- function(med_bf_object,
                            annot) {
  
  bf_dat <- data.frame(protein.id = names(med_bf_object$lnBF_med), lnBF = med_bf_object$lnBF_med) %>%
    dplyr::left_join(annot)
  
  approx_pval <- nrow(bf_dat)
  for (i in 1:nrow(bf_dat)) {
    approx_pval[i] <- mean(c(bf_dat$lnBF_med[i] < (bf_dat %>% 
                                                     dplyr::filter(chr != bf_dat$chr[i]) %>% 
                                                     dplyr::pull(lnBF_med)), TRUE))
  }
  bf_dat$approx_pval <- approx_pval
  bf_dat
}

#' Function that calculate permutation-based p-values for Bayes factors
#'
#' This function takes Bayes factor results from mediation_bf() and calculates permutation p-values based
#' on permuting the mediator.
#'
#' @export
#' @examples get_perm_pval()
get_perm_pval <- function(y, 
                          M, 
                          X, 
                          Z = NULL,
                          num_perm = 1000,
                          verbose = FALSE,
                          bf_type = "lnBF_totmed",
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
    sub_Z <- Z[apply(Z, 1, function(x) paste(x, collapse = "")) == paste(within_cov[i,], collapse = ""),, drop = FALSE]
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
                               verbose = verbose,
                               ...)
      perm_bf_mat[i, j] <- perm_med[[bf_type]]
    } 
    print(paste("Perm", i, "done out of", num_perm))
  }
  perm_pval <- sapply(1:ncol(perm_bf_mat), function(i) mean(actual_bf[[bf_type]][i] < perm_bf_mat[,i]))
  bf_dat <- data.frame(mediator = colnames(M), 
                       bf = actual_bf[[bf_type]], 
                       perm_pval = perm_pval, 
                       thresh = apply(perm_bf_mat, 2, function(x) quantile(x, probs = 0.95)))
  names(bf_dat)[which(names(bf_dat) == "bf")] <- bf_type
  bf_dat
}
