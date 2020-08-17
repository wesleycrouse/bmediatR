#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`

# Conversion function taken from qtl2convert
map_df_to_list <- function (map, chr_column = "chr", pos_column = "cM", marker_column = "marker", 
          Xchr = c("x", "X")) {
  
  if (is.null(marker_column)) {
    marker_column <- "qtl2tmp_marker"
    map[, marker_column] <- rownames(map)
  }
  if (!(marker_column %in% colnames(map))) 
    stop("Column \"", marker_column, "\" not found.")
  if (!(chr_column %in% colnames(map))) 
    stop("Column \"", chr_column, "\" not found.")
  if (!(pos_column %in% colnames(map))) 
    stop("Column \"", pos_column, "\" not found.")
  marker <- map[, marker_column]
  chr <- map[, chr_column]
  uchr <- unique(chr)
  pos <- map[, pos_column]
  result <- split(as.numeric(pos), factor(chr, levels = uchr))
  marker <- split(marker, factor(chr, levels = uchr))
  for (i in seq(along = result)) names(result[[i]]) <- marker[[i]]
  is_x_chr <- rep(FALSE, length(result))
  names(is_x_chr) <- names(result)
  if (!is.null(Xchr)) {
    Xchr_used <- Xchr %in% names(is_x_chr)
    if (any(Xchr_used)) {
      Xchr <- Xchr[Xchr_used]
      is_x_chr[Xchr] <- TRUE
    }
  }
  attr(result, "is_x_chr") <- is_x_chr
  result
}

#' Bayes factor genome plot function
#'
#' This function takes Bayes factor results from mediation_bf() and plots the genome-wide mediation scan.
#'
#' @param med_bf_object Output from mediation_bf(). 
#' @param bf_type DEFAULT: "lnBF_med". Bayes factor to be displayed. 
#' @param med_annot Annotation data for -omic mediators.
#' @param include_chr DEFAULT: c(1:19, "X"). Chromosomes to include in plot.
#' @param expland_lim_factor DEFAULT: 0.025. Scale to increase plot limits by.
#' @param label_thresh DEFAULT: NULL. Label mediators that surpass label_thresh. Default does not add labels.
#' @param qtl_dat DEFAULT: NULL. QTL data that includes position of QTL and outcome. Adds ticks to the figure.
#' @export
#' @examples plot_bf()
plot_bf <- function(med_bf_object, 
                    bf_type = "lnBF_med",
                    med_annot, 
                    include_chr = c(1:19, "X"), 
                    expand_lim_factor = 0.025, 
                    label_thresh = NULL, 
                    bgcol = "white", altcol = "gray", altbgcol = "white", hlines_col = "gray80", col = "black", cex = 0.75,
                    qtl_dat = NULL,
                    outcome_symbol = NULL,
                    ...) {
  
  bf_type <- bf_type[1]
  
  bf <- matrix(med_bf_object[[bf_type]], ncol = 1)
  rownames(bf) <- names(med_bf_object[[bf_type]])
  class(bf) <- "scan1"
  
  med_map_df <- med_annot %>%
    dplyr::select(protein.id, symbol, chr, middle) %>%
    dplyr::filter(chr %in% include_chr) %>%
    dplyr::mutate(chr = factor(chr, levels = c(1:19, "X"))) %>%
    as.data.frame %>% 
    dplyr::arrange(chr)
  if (!is.null(qtl_dat)) {
    ## Add QTL to map for plotting
    med_map_df <- dplyr::bind_rows(med_map_df,
                                   data.frame(protein.id = "QTL", symbol = "QTL", chr = qtl_dat$chr, middle = qtl_dat$pos))
    
  }
  med_map <- map_df_to_list(map = med_map_df, marker_column = "protein.id", pos_column = "middle")
  
  gap <- sum(qtl2::chr_lengths(map))/100
  
  lim_shift <- (max(bf[,1]) - min(bf[,1])) * expand_lim_factor
  qtl2:::plot.scan1(bf, map = med_map, ylab = "Log Bayes factor", type = "p", pch = 20, ylim = c(min(bf[,1]) - lim_shift, max(bf[,1]) + lim_shift),
       bgcol = bgcol, altcol = altcol, altbgcol = altbgcol, hlines_col = hlines_col, col = col, cex = cex, gap = gap,
       ...)
  
  xpos <- qtl2:::map_to_xpos(map = med_map, gap = gap)
  
  if (!is.null(label_thresh) & any(bf > label_thresh)) {
    labels <- rownames(bf)[bf > label_thresh]
    
    label_map_df <- med_map_df %>%
      filter(protein.id %in% labels) 
    
    for (i in 1:nrow(label_map_df)) {
      lab_pos <- xpos[label_map_df$protein.id[i]]
      lab_bf <- bf[label_map_df$protein.id[i],]
      
      text(x = lab_pos, y = lab_bf, label_map_df$symbol[i], font = 3)
    }
  }
  if (!is.null(outcome_symbol)) {
    rug(x = xpos[med_annot %>% 
                   dplyr::filter(symbol == outcome_symbol) %>% pull(protein.id)],
        lwd = 3,
        col = "black")
  }
  if (!is.null(qtl_dat)) {
    rug(x = xpos["QTL"],
        lwd = 3,
        col = "red")
  }
}

#' Posterior probability genome plot function
#'
#' This function takes posterior probability results from mediation_bf() and plots the genome-wide mediation scan.
#'
#' @export
#' @examples plot_posterior()
plot_posterior <- function(med_bf_object, 
                           post_col = c(4, 8),
                           med_annot, 
                           include_chr = c(1:19, "X"), 
                           expand_lim_factor = 0.025, 
                           label_thresh = NULL, 
                           bgcol = "white", altcol = "gray", altbgcol = "white", hlines_col = "gray80", col = "black", cex = 0.75,
                           qtl_dat = NULL,
                           outcome_symbol = NULL,
                           ...) {
  
  post_p <- matrix(apply(exp(med_bf_object$ln_post_c[,c(4,8)]), 1, function(x) sum(x)), ncol = 1)
  rownames(post_p) <- rownames(med_bf_object$ln_post_c)
  class(post_p) <- "scan1"
  
  med_map_df <- med_annot %>%
    dplyr::select(protein.id, symbol, chr, middle) %>%
    dplyr::filter(chr %in% include_chr) %>%
    dplyr::mutate(chr = factor(chr, levels = c(1:19, "X"))) %>%
    as.data.frame %>% 
    dplyr::arrange(chr)
  if (!is.null(qtl_dat)) {
    ## Add QTL to map for plotting
    med_map_df <- dplyr::bind_rows(med_map_df,
                                   data.frame(protein.id = "QTL", symbol = "QTL", chr = qtl_dat$chr, middle = qtl_dat$pos))
    
  }
  med_map <- map_df_to_list(map = med_map_df, marker_column = "protein.id", pos_column = "middle")
  
  gap <- sum(qtl2::chr_lengths(map))/100
  
  lim_shift <- (max(post_p[,1]) - min(post_p[,1])) * expand_lim_factor
  qtl2:::plot.scan1(post_p, map = med_map, ylab = "Posterior Probability", type = "p", pch = 20, ylim = c(min(post_p[,1]) - lim_shift, max(post_p) + lim_shift),
                    bgcol = bgcol, altcol = altcol, altbgcol = altbgcol, hlines_col = hlines_col, col = col, cex = cex, gap = gap,
                    ...)
  
  xpos <- qtl2:::map_to_xpos(map = med_map, gap = gap)
  
  if (!is.null(label_thresh) & any(post_p > label_thresh)) {
    labels <- rownames(post_p)[post_p > label_thresh]
    
    label_map_df <- med_map_df %>%
      filter(protein.id %in% labels) 
    
    for (i in 1:nrow(label_map_df)) {
      lab_pos <- xpos[label_map_df$protein.id[i]]
      lab_post_p <- post_p[label_map_df$protein.id[i],]
      
      text(x = lab_pos, y = lab_post_p, label_map_df$symbol[i], font = 3)
    }
  }
  if (!is.null(outcome_symbol)) {
    rug(x = xpos[med_annot %>% 
                   dplyr::filter(symbol == outcome_symbol) %>% pull(protein.id)],
        lwd = 3,
        col = "black")
  }
  if (!is.null(qtl_dat)) {
    rug(x = xpos["QTL"],
        lwd = 3,
        col = "red")
  }
}

