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

#' Barplot of posterior model probabilities function
#'
#' This function takes posterior probability results from mediation_bf() and plots barplots of posterior model probabilities.
#'
#' @param bmediatR_object Output from bmediatR(). 
#' @param med_annot Annotation data for -omic mediators.
#' @param mediator_id Which mediator to plot posterior model probabilities for.
#' @param med_var DEFAULT: "protein.id". The column in med_annot to be used as a mediator id.
#' @param stack DEFAULT: FALSE. If TRUE, a stacked barplot is produced. If FALSE, a staggered
#' barplot is produced.
#'
#' @export
#' @examples plot_posterior_bar()
plot_posterior_bar <- function(bmediatR_object,
                               med_annot = NULL,
                               mediator_id,
                               med_var = "protein.id",
                               stack = FALSE,
                               bar_col = c("seagreen4", "seagreen1", "skyblue", "goldenrod1", "goldenrod4", "gray"),
                               relabel_x = NULL,
                               main = NULL) {
  
  ## Flag for reactive model
  prior_mat <- bmediatR_object$ln_prior_c
  model_flag <- is.finite(prior_mat)
  names(model_flag) <- colnames(prior_mat)
  
  ## Long names
  long_names <- c("other non-med", "other non-med", "other non-med", 
                  "complete med", "other non-med", "other non-med",
                  "co-local", "partial med", "other non-med",
                  "complete med (react)", "other non-med", "partial med (react)")
  names(long_names) <- colnames(prior_mat)
  

  bar_col <- bar_col[c(model_flag[c("0,1,1", "1,1,1", "1,1,0", "1,1,*", "1,0,*")], TRUE)]

  posterior_dat <- exp(bmediatR_object$ln_post_c) %>%
    as.data.frame %>%
    tibble::rownames_to_column(med_var) %>%
    dplyr::rename("partial med" = `1,1,1`,
                  "complete med" = `0,1,1`,
                  "co-local" = `1,1,0`,
                  "partial med (react)" = `1,1,*`,
                  "complete med (react)" = `1,0,*`)

  # Using annotation information
  if (!is.null(med_annot)) {
    posterior_dat <- posterior_dat %>%
      dplyr::left_join(med_annot %>%
                         dplyr::select(tidyselect::all_of(med_var), symbol))
  } else {
    posterior_dat <- posterior_dat %>%
      dplyr::mutate(symbol = get(med_var))
  }
  
  ## Calculating non-mediation or co-local probability
  posterior_dat <- posterior_dat %>%
    dplyr::left_join(posterior_dat %>%
                       dplyr::select(tidyselect::all_of(med_var), contains(",")) %>%
                       dplyr::mutate("other non-med" = rowSums(.[-1]))) %>%
    dplyr::select(-contains(",")) %>%
    tidyr::gather(key = model, value = post_p, -c(tidyselect::all_of(med_var), symbol))

  ## Set factors
  models_use <- unique(long_names[model_flag])
  posterior_dat <- posterior_dat %>%
    dplyr::filter(model %in% models_use) %>%
    dplyr::mutate(model = factor(model, levels = c("complete med", "partial med", "co-local", "partial med (react)", "complete med (react)", "other non-med")))
  
  bar_theme <- ggplot2::theme(panel.grid.major = ggplot2::element_blank(), 
                              panel.grid.minor = ggplot2::element_blank(),
                              panel.background = ggplot2::element_blank(), 
                              axis.line = ggplot2::element_line(colour = "black"),
                              plot.title = ggplot2::element_text(hjust = 0.5, size = 16, face = "plain"), 
                              axis.title.x = ggplot2::element_blank(),
                              axis.text.x = ggplot2::element_text(hjust = 0.5, size = 14, face = "plain"),
                              axis.title.y = ggplot2::element_text(size = 14, face = "plain"),
                              axis.text.y = ggplot2::element_text(size = 14, face = "plain"),
                              axis.ticks.y = ggplot2::element_blank(),
                              legend.title = ggplot2::element_text(size = 14),
                              legend.text = ggplot2::element_text(size = 14))
  
  if (!is.null(relabel_x)) {
    posterior_dat$symbol <- relabel_x
  }
  
  p <- ggplot2::ggplot(data = posterior_dat %>% 
                         dplyr::filter((!!as.symbol(med_var)) == mediator_id)) +
    ggplot2::scale_fill_manual(values = bar_col) +
    ggplot2::ylab("Posterior model probability") + 
    bar_theme
  if (!is.null(main)) {
    p <- p + ggplot2::ggtitle(main)
  }
  if (stack) {
    p <- p + ggplot2::geom_col(ggplot2::aes(x = symbol, y = post_p, fill = model), width = 0.5) 
  } else {
    p <- p + ggplot2::geom_bar(ggplot2::aes(x = symbol, y = post_p, fill = model), position = "dodge", stat = "summary", fun = "mean") +
      ggplot2::geom_hline(yintercept = c(0, 1), col = "gray", linetype = "dashed")
  }
  
  p
}

#' Posterior odds genome plot function
#'
#' This function takes the posterior odds results from bmediatR() and plots the genome-wide scan.
#'
#' @param bmediatR_object Output from bmediatR(). 
#' @param model_type DEFAULT: "mediation". Specifies which model(s)'s posterior probabilities are to be included in the numerator of the posterior odds and then displayed for
#' for genome-wide mediators. 
#' @param med_annot Annotation data for -omic mediators.
#' @param include_chr DEFAULT: c(1:19, "X"). Chromosomes to include in plot.
#' @param expland_lim_factor DEFAULT: 0.025. Scale to increase plot limits by.
#' @param label_thresh DEFAULT: NULL. Label mediators that surpass label_thresh. Default does not add labels.
#' @param label_only_chr DEFAULT: NULL. Only label mediators that pass label_thresh on the specified chromosome.
#' @param qtl_dat DEFAULT: NULL. QTL data that includes position of QTL and outcome. Adds ticks to the figure.
#' @export
#' @examples plot_posterior_odds()
plot_posterior_odds <- function(bmediatR_object, 
                                model_type = c("mediation", "partial", "complete", "colocal"),
                                med_annot, 
                                med_var = "protein.id",
                                include_chr = c(1:19, "X"), 
                                expand_lim_factor = 0.025, 
                                label_thresh = NULL, 
                                label_only_chr = NULL,
                                bgcol = "white", altcol = "gray", altbgcol = "white", 
                                hlines_col = "gray80", col = "black", cex = 0.75,
                                qtl_dat = NULL,
                                outcome_symbol = NULL,
                                ...) {
  
  model_type <- model_type[1]
  
  post_odds <- matrix(bmediatR_object[["ln_post_odds"]][,model_type], ncol = 1)
  rownames(post_odds) <- rownames(bmediatR_object[["ln_post_odds"]])
  class(post_odds) <- "scan1"
  
  med_map_df <- med_annot %>%
    dplyr::select(tidyselect::all_of(med_var), symbol, chr, middle) %>%
    dplyr::filter(chr %in% include_chr) %>%
    dplyr::mutate(chr = factor(chr, levels = c(1:19, "X"))) %>%
    as.data.frame %>% 
    dplyr::arrange(chr)
  if (!is.null(qtl_dat)) {
    ## Add QTL to map for plotting
    med_map_df <- dplyr::bind_rows(med_map_df,
                                   qtl_dat %>%
                                     dplyr::mutate((!!as.symbol(med_var)) := "QTL",
                                                   symbol = "QTL") %>%
                                     dplyr::rename(middle = pos) %>%
                                     dplyr::select(tidyselect::all_of(med_var), symbol, chr, middle))
  }
  med_map <- map_df_to_list(map = med_map_df, marker_column = med_var, pos_column = "middle")
  
  gap <- sum(qtl2::chr_lengths(map))/100
  
  lim_shift <- (max(post_odds[,1]) - min(post_odds[,1])) * expand_lim_factor
  qtl2:::plot.scan1(post_odds, map = med_map, ylab = "Log posterior odds", type = "p", pch = 20, 
                    ylim = c(min(post_odds[,1]) - lim_shift, max(post_odds[,1]) + lim_shift),
                    bgcol = bgcol, altcol = altcol, altbgcol = altbgcol, hlines_col = hlines_col, col = col, 
                    cex = cex, gap = gap,
                    ...)
  
  xpos <- qtl2:::map_to_xpos(map = med_map, gap = gap)
  
  ## Mediator labels
  label_dat <- matrix(bmediatR_object[["ln_post_odds"]][,model_type], ncol = 1)
  colnames(label_dat) <- "post_odds"
  rownames(label_dat) <- rownames(bmediatR_object[["ln_post_odds"]])
  label_dat <- label_dat %>%
    as.data.frame %>%
    tibble::rownames_to_column(med_var) %>%
    dplyr::left_join(med_map_df)
  if (!is.null(label_only_chr)) {
    label_dat <- label_dat %>%
      dplyr::filter(chr == label_only_chr)
  } else {
    label_dat <- label_dat %>%
      dplyr::filter(chr %in% include_chr)
  }
  label_post_odds <- label_dat %>%
    dplyr::select(tidyselect::all_of(med_var), post_odds) %>%
    tibble::column_to_rownames(med_var) %>%
    as.matrix()
  
  if (!is.null(label_thresh) & any(label_post_odds > label_thresh)) {
    labels <- rownames(label_post_odds)[label_post_odds > label_thresh]
    
    label_map_df <- med_map_df %>%
      filter((!!as.symbol(med_var)) %in% labels) 
    
    for (i in 1:nrow(label_map_df)) {
      lab_pos <- xpos[label_map_df[i, med_var]]
      lab_post_odds <- post_odds[label_map_df[i, med_var],]
      
      text(x = lab_pos, y = lab_post_odds, label_map_df$symbol[i], font = 3)
    }
  }
  if (!is.null(outcome_symbol)) {
    rug(x = xpos[med_annot %>% 
                   dplyr::filter(symbol == outcome_symbol) %>% 
                   pull(tidyselect::all_of(med_var))],
        lwd = 3,
        col = "black")
  }
  if (!is.null(qtl_dat)) {
    rug(x = xpos["QTL"],
        lwd = 3,
        col = "red")
  }
}


