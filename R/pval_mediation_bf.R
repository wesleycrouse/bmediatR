#' Function that calculate approximate p-values for posterior odds
#'
#' This function takes posterior odds results from bmediatR() and calculates approximate p-values based
#' on comparisons to genome results, using the assumption that the vast majority of candidates are not
#' mediators and can thus be used to characterize a null distribution of posterior odds.
#'
#' @param bmediatR_object Output from bmediatR(). 
#' @param model_type DEFAULT: "mediation". Specifies which model(s)'s posterior probabilities are to be included in the numerator of the posterior odds and then displayed for
#' for genome-wide mediators. 
#' @param med_annot Annotation data for -omic mediators.
#'
#' @export
#' @examples get_approx_pval()
get_approx_pval <- function(bmediatR_object,
                            model_type = c("mediation", "partial", "complete", "colocal"),
                            med_annot,
                            med_var = "protein.id") {
  
  model_type <- model_type[1]
  
  post_odds_dat <- bmediatR_object$ln_post_odds %>%
    as.data.frame %>%
    dplyr::select(tidyselect::all_of(model_type)) %>%
    tibble::rownames_to_column(med_var) %>%
    dplyr::left_join(med_annot)
  
  approx_pval <- nrow(post_odds_dat)
  for (i in 1:nrow(post_odds_dat)) {
    approx_pval[i] <- mean(c(post_odds_dat[i,model_type] < (post_odds_dat %>% 
                                                              dplyr::filter(chr != post_odds_dat$chr[i]) %>% 
                                                              dplyr::pull(tidyselect::all_of(model_type))), TRUE))
  }
  post_odds_dat$approx_pval <- approx_pval
  post_odds_dat <- post_odds_dat %>%
    dplyr::rename(ln_post_odds = tidyselect::all_of(model_type))
  
  post_odds_dat
}

#' Function that calculate permutation-based p-values for posterior odds
#'
#' This function takes posterior odds results from bmediatR() and calculates permutation p-values based
#' on permuting the mediator.
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
#'
#' @export
#' @examples get_perm_pval()
get_perm_pval <- function(y, M, X, Z = NULL, w = NULL,
                          model_type = c("mediation", "partial", "complete", "colocal"),
                          ln_prior_c = "complete",
                          num_perm = 1000,
                          verbose = FALSE,
                          ...) {
  
  actual_po <- bmediatR(y = y, M = M, X = X, Z = Z, w = w,
                        ln_prior_c = ln_prior_c,
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
  perm_po_mat <- matrix(NA, nrow = num_perm, ncol = ncol(M))
  colnames(perm_po_mat) <- colnames(M)
  for(i in 1:num_perm) {
    ## Rename things
    perm_M <- M
    rownames(perm_M) <- as.character(perm_mat[,i])
    for (j in 1:ncol(M)) {
      perm_med <- bmediatR(y = y, M = perm_M[,j,drop = FALSE][rownames(M),], X = X, Z = Z, w = w,
                           ln_prior_c = ln_prior_c,
                           verbose = verbose,
                           ...)
      perm_po_mat[i, j] <- perm_med$ln_post_odds[,model_type]
    } 
    print(paste("Perm", i, "done out of", num_perm))
  }
  perm_pval <- sapply(1:ncol(perm_po_mat), function(i) mean(actual_po$ln_post_odds[i, model_type] < perm_po_mat[,i]))
  post_odds_dat <- data.frame(mediator = colnames(M), 
                              ln_post_odds = actual_po$ln_post_odds[,model_type], 
                              perm_pval = perm_pval, 
                              thresh = apply(perm_po_mat, 2, function(x) quantile(x, probs = 0.95)))
  post_odds_dat
}
