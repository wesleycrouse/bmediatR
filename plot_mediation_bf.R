library(tidyverse)

## Plot function
plot_bf <- function(med_bf_object, 
                    bf_type = c("lnBF_med", "lnBF_coloc"),
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
  rownames(bf) <- names(med_bf_object[bf_type])
  class(bf) <- "scan1"
  
  med_map_df <- med_annot %>%
    dplyr::select(protein.id, symbol, chr, middle) %>%
    filter(chr %in% include_chr) %>%
    mutate(chr = factor(chr, levels = c(1:19, "X"))) %>%
    as.data.frame %>% 
    arrange(chr)
  if (!is.null(qtl_dat)) {
    ## Add QTL to map for plotting
    med_map_df <- bind_rows(med_map_df,
                            data.frame(protein.id = "QTL", symbol = "QTL", chr = qtl_dat$chr, middle = qtl_dat$pos))
    
  }
  med_map <- qtl2convert::map_df_to_list(map = med_map_df, marker_column = "protein.id", pos_column = "middle")
  
  gap <- sum(chr_lengths(map))/100
  
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
    rug(x = xpos[med_annot %>% filter(symbol == outcome_symbol) %>% pull(protein.id)],
        lwd = 3,
        col = "black")
  }
  if (!is.null(qtl_dat)) {
    rug(x = xpos["QTL"],
        lwd = 3,
        col = "red")
  }
}
