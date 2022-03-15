#################################################################################################################
#################################################################################################################
######
######                Title: bmediatr_manuscript_support_functions.R
######                            
######                Description: Analysis and plotting functions used to generate results reported in 
######                             the bmediatR manuscript
######     
######                            
######                Manuscript: Crouse & Keele et al. 
######
######                Author: Greg Keele
######
#################################################################################################################
#################################################################################################################

########################################################
###### 
###### Nice scan plot
######
########################################################
# Pulled from qtl2convert
map_list_to_df <- function(map_list, chr_column = "chr", pos_column = "pos", marker_column = "marker") {
  
  nmar <- vapply(map_list, length, 1)
  markers <- unlist(lapply(map_list, names))
  result <- data.frame(chr = rep(names(map_list), nmar), pos = unlist(map_list), 
                       marker = markers, stringsAsFactors = FALSE)
  rownames(result) <- markers
  names(result)[1] <- chr_column
  names(result)[2] <- pos_column
  if (is.null(marker_column)) 
    result <- result[, -3, drop = FALSE]
  else names(result)[3] <- marker_column
  result
}

nice_scan <- function(qtl_scan,
                      med_scan = NULL,
                      med_highlight_dat = NULL,
                      map,
                      bgcolor = "white",
                      altbgcolor = "white",
                      bgcolor_alpha = 0.6,
                      altbgcolor_alpha = 0.6,
                      color,
                      altcolor,
                      medcol = "gray90",
                      outlinecol = "black",
                      hline = NULL,
                      annot_dat,
                      target_symbol = NULL,
                      chr = c(1:19, "X"),
                      main = "",
                      my_y_max = NULL,
                      rug_col = "black") {
  
  include_chr <- chr
  if (is.null(med_scan)) {
    y_max <- max(qtl_scan)
  } else {
    med_scan <- med_scan %>% 
      filter(chr %in% chr)
    y_max <- max(c(qtl_scan, med_scan$lod))
  }
  if (!is.null(my_y_max)) {
    y_max <- max(y_max, my_y_max)
  }
  plot(qtl_scan, map, 
       ylim = c(0, y_max), 
       bgcolor = bgcolor, 
       altbgcolor = altbgcolor, 
       col = "white", 
       gap = 0,
       chr = chr,
       main = main)
  map_df <- map_list_to_df(map)
  j <- 1
  for (i in chr) {
    x <- xpos_scan1(map = map, thepos = map_df %>% filter(chr == i) %>% pull(pos), thechr = map_df %>% filter(chr == i) %>% pull(chr), gap = 0, chr = chr)
    y <- qtl_scan[map_df %>% filter(chr == i) %>% pull(marker),]
    polygon(x = c(x[1], x, x[length(x)]),
            y = c(0, y, 0), 
            border = NA,
            col = ifelse(j %% 2 != 0, 
                         scales::alpha(color, bgcolor_alpha), 
                         scales::alpha(altcolor, altbgcolor_alpha)))
    j <- j + 1
  }
  if (!is.null(med_scan)) {
    points(x = qtl2::xpos_scan1(map, 
                                thechr = med_scan$chr, 
                                thepos = med_scan$pos, 
                                gap = 0,
                                chr = chr), 
           y = med_scan$lod,
           col = medcol, 
           pch = 20)
  }
  
  plot(qtl_scan, map, 
       ylim = c(0, y_max), 
       bgcolor = bgcolor, 
       altbgcolor = altbgcolor, 
       col = outlinecol, 
       lwd = 0.5,
       gap = 0,
       add = TRUE,
       chr = chr)
  if (!is.null(med_highlight_dat)) {
    for (i in 1:nrow(med_highlight_dat)) {
      points(x = qtl2::xpos_scan1(map, 
                                  thechr = med_scan %>% filter(symbol == med_highlight_dat$symbol[i]) %>% pull(chr), 
                                  thepos = med_scan %>% filter(symbol == med_highlight_dat$symbol[i]) %>% pull(pos), 
                                  gap = 0,
                                  chr = chr), 
             y = med_scan %>% filter(symbol == med_highlight_dat$symbol[i]) %>% pull(lod),
             col = med_highlight_dat$col[i], 
             bg = med_highlight_dat$bgcol[i],
             cex = med_highlight_dat$cex[i],
             pch = med_highlight_dat$pch[i])
    }
  }
  if (!is.null(target_symbol)) {
    for (i in 1:length(target_symbol)) {
      rug(x = qtl2::xpos_scan1(map, 
                               thechr = annot_dat %>% filter(symbol == target_symbol[i]) %>% pull(chr),
                               thepos = annot_dat %>% filter(symbol == target_symbol[i]) %>% pull(middle),
                               gap = 0,
                               chr = chr),
          col = rug_col,
          lwd = 3)
    }
  }
  if (!is.null(hline)) {
    abline(h = hline, lty = 2)
  }
}

######################################################################
###### 
###### Function to calculate BLUP allele effect for lodpeaks table
######
######################################################################
calc_allele_effects_qtl2 <- function(qtl_table,
                                     genoprobs,
                                     phenotype_mat,
                                     covar_mat = NULL,
                                     K = NULL,
                                     data.col = c("protein.id", "gene.id"),
                                     map,
                                     add.marker.id = TRUE) {
  
  data.col <- data.col[1]
  effect_mat <- effect_se_mat <- matrix(NA, nrow = nrow(qtl_table), ncol = 8)
  colnames(effect_mat) <- LETTERS[1:8]
  colnames(effect_se_mat) <- paste(LETTERS[1:8], "SE", sep = "_")
  if (add.marker.id) {
    qtl_markers <- rep(NA, nrow(qtl_table))
  }
  for (i in 1:nrow(qtl_table)) {
    this_map <- map[[as.character(qtl_table$chr[i])]]
    
    if (add.marker.id) {
      this_marker <- find_marker(map = map, pos = qtl_table$pos[i], chr = qtl_table$chr[i])
      qtl_markers[i] <- this_marker
    } else {
      this_marker <- qtl_table$marker.id[i]
    }
    
    this_genoprobs <- qtl2::pull_genoprobpos(genoprobs, 
                                             marker = this_marker)
    
    this_fit <- qtl2::fit1(genoprobs = this_genoprobs, 
                           pheno = phenotype_mat[,qtl_table[i, data.col], drop = FALSE], 
                           addcovar = covar_mat,
                           kinship = K[[as.character(qtl_table$chr[i])]],
                           blup = TRUE)
    if (all(this_fit$coef[LETTERS[1:8]] == 0)) {
      this_fit <- qtl2::fit1(genoprobs = this_genoprobs, 
                             pheno = phenotype_mat[,qtl_table[i, data.col], drop = FALSE], 
                             addcovar = covar_mat,
                             kinship = K[[as.character(qtl_table$chr[i])]],
                             blup = FALSE)
    }
    effect_mat[i,] <- this_fit$coef[LETTERS[1:8]]
    effect_se_mat[i,] <- this_fit$SE[LETTERS[1:8]]
  }
  if (add.marker.id) {
    updated_qtl_table <- bind_cols(bind_cols(qtl_table, 
                                             data.frame(marker.id = qtl_markers)), 
                                   as.data.frame(effect_mat),
                                   as.data.frame(effect_se_mat))
  } else {
    updated_qtl_table <- bind_cols(qtl_table, 
                                   as.data.frame(effect_mat),
                                   as.data.frame(effect_se_mat))
  }
}

######################################################################
###### 
###### Plot comparing haplotype effects at QTL with error bars
######
######################################################################
effects_w_ci_plot <- function(effect_dat,
                              use_se = TRUE,
                              slope_multiplier = 1,
                              col = as.character(qtl2::CCcolors),
                              xlab = "QTL effects on Y",
                              ylab = "QTL effects on M",
                              main = "") {
  
  gg_theme <- theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(), 
                    axis.line = element_line(colour = "black"),
                    plot.title = element_text(hjust = 0.5), 
                    axis.text = element_text(size = 12),
                    axis.title = element_text(size = 12))
  
  multiplier <- ifelse(use_se, 1, 1.96)
  xmin = min(effect_dat$y_effect - multiplier*effect_dat$y_SE)
  xmax = max(effect_dat$y_effect + multiplier*effect_dat$y_SE)
  ymin = min(effect_dat$m_effect - multiplier*effect_dat$m_SE)
  ymax = max(effect_dat$m_effect + multiplier*effect_dat$m_SE)
  
  p <- ggplot(data = effect_dat, aes(x = y_effect, y = m_effect, col = strain)) + 
    geom_point() +
    geom_errorbar(aes(ymin = m_effect - multiplier*m_SE, ymax = m_effect + multiplier*m_SE), width = 0, size = 1.001, alpha = 0.7) + 
    geom_errorbarh(aes(xmin = y_effect - multiplier*y_SE, xmax = y_effect + multiplier*y_SE), height = 0, size = 1.001, alpha = 0.7) + 
    scale_color_manual(values = col) + 
    xlim(xmin, xmax) + ylim(ymin, ymax) +
    geom_vline(xintercept = 0, col = "gray", size = 0.5, linetype = "dotted") + 
    geom_hline(yintercept = 0, col = "gray", size = 0.5, linetype = "dotted") + 
    geom_abline(slope = 1*slope_multiplier, intercept = 0, linetype = "dashed", size = 0.5) + 
    xlab(xlab) + ylab(ylab) +
    annotate("text",
             label = paste("r =", round(cor(effect_dat$m_effect, effect_dat$y_effect), 3)), 
             size = 4,
             x = xmin + (xmax - xmin)/2, 
             y = ymax) +
    ggtitle(main) +
    gg_theme +
    guides(col = "none")
  p
}

######################################################################
###### 
###### Calculate the proportion of variance of y explained by X
######
######################################################################
calc_effect_size <- function(y, X,
                             Z = NULL,
                             w = NULL,
                             make_intercept = TRUE) {
  
  if (is.matrix(y)) {
    y <- y[,1]
  }
  
  if (make_intercept) {
    int_mat <- matrix(rep(1, length(y)), ncol = 1); rownames(int_mat) <- names(y)
  } else{
    int_mat <- NULL
  }
  
  fit_alt <- lm.wfit(x = cbind(int_mat, Z, X), y = y, w = w)
  fit_null <- lm.wfit(x = cbind(int_mat, Z), y = y, w = w)
  
  #effect_size <- 1 - sum((y[names(fit_alt$fitted.values)] - fit_alt$fitted.values)^2)/sum((y[names(fit_alt$fitted.values)] - fit_null$fitted.values)^2)
  effect_size <- 1 - sum((y - fit_alt$fitted.values)^2)/sum((y - fit_null$fitted.values)^2)
  
  effect_size
}

######################################################################
###### 
###### Calculate empirical hyper priors for 
###### X -> M, X -> Y, and M -> Y
######
######################################################################
calc_trio_effect_sizes <- function(y, m, X,
                                   Z = NULL, Z_y = NULL, Z_m = NULL,
                                   w = NULL, w_y = NULL, w_m = NULL,
                                   make_intercept = TRUE,
                                   align_data = TRUE) {
  
  processed_data <- bmediatR:::process_data(y = y, M = m, X = X,
                                            Z_y = Z_y, Z_M = Z_m,
                                            w_y = w_y, w_M = w_m, 
                                            align_data = align_data,
                                            verbose = FALSE)
  
  y <- processed_data$y
  m <- processed_data$M[,1]
  X <- processed_data$X
  Z_y <- processed_data$Z_y
  Z_m <- processed_data$Z_M
  w_y <- processed_data$w_y
  w_m <- processed_data$w_M
  
  x_m <- calc_effect_size(Z = Z_m, 
                          y = m, 
                          X = X, 
                          w = w_m, 
                          make_intercept = make_intercept)
  m_y <- calc_effect_size(Z = Z_y, 
                          y = y, 
                          X = m, 
                          w = w_y, 
                          make_intercept = make_intercept)
  x_y <- calc_effect_size(Z = Z_y,
                          y = y, 
                          X = X, 
                          w = w_y, 
                          make_intercept = make_intercept)
  
  results <- c(x_m, m_y, x_y)
  names(results) <- c("x_m", "m_y", "x_y")
  results
}

renormalize_effect_size_ratio <- function (effect_size,
                                           min_noise = 0.001) {
  
  effect_size <- effect_size/(1 + min_noise)
  ratio <- effect_size/(1 - effect_size)
  ratio
}

######################################################################
###### 
###### Function to run local scan
######
######################################################################
run_local_scan <- function(annot,
                           genoprobs,
                           phenotype_mat,
                           covar_mat = NULL,
                           K = NULL,
                           data.col = c("protein.id", "gene.id"),
                           window = 0,
                           map,
                           markers) {
  
  data.col <- data.col[1]
  lod_vec <- marker_vec <- rep(NA, nrow(annot))
  names(lod_vec) <- annot[,data.col, drop = TRUE]
  names(marker_vec) <- annot[,data.col, drop = TRUE]

  for (i in 1:nrow(annot)) {
    if (!annot$chr[i] %in% c("MT", "Y")) {
      if (window == 0) {
        marker_vec[i] <- markers %>%
          filter(chr == annot$chr[i]) %>%
          mutate(dist = abs(pos - (annot$start[i] + (annot$end[i] - annot$start[i])/2))) %>%
          filter(dist == min(dist)) %>%
          pull(marker.id)
        this_genoprobs <- qtl2::pull_genoprobpos(genoprobs, 
                                                 marker = marker_vec[i])
        
        this_fit <- qtl2::fit1(genoprobs = this_genoprobs, 
                               pheno = phenotype_mat[,annot[i, data.col, drop = TRUE], drop = FALSE], 
                               addcovar = covar_mat,
                               kinship = K[[as.character(annot$chr[i])]],
                               blup = FALSE)
        
        lod_vec[i] <- this_fit$lod
      } else {
        these_genoprobs <- qtl2::pull_genoprobint(genoprobs = genoprobs, 
                                                  map = map, 
                                                  chr = annot$chr[i],
                                                  interval = c(annot$start[i] + (annot$end[i] - annot$start[i])/2 - window,
                                                               annot$start[i] + (annot$end[i] - annot$start[i])/2 + window))
        these_fits <- qtl2::scan1(genoprobs = these_genoprobs, 
                                  pheno = phenotype_mat[,annot[i, data.col, drop = TRUE], drop = FALSE], 
                                  addcovar = covar_mat,
                                  kinship = K[[as.character(annot$chr[i])]])
        this_max <- which.max(these_fits)
        lod_vec[i] <- these_fits[this_max]
        marker_vec[i] <- rownames(these_fits)[this_max]
      }
      print(i)
    }
  }
  
  local_lod_table <- data.frame(holder = names(lod_vec), lod = lod_vec, marker.id = marker_vec)
  names(local_lod_table)[1] <- data.col
    
  local_lod_table <- local_lod_table %>% 
    left_join(annot %>%
                dplyr::select(tidyselect::all_of(data.col), chr, start, end) %>%
                dplyr::rename(gene.start = start,
                              gene.end = end) %>%
                mutate(gene.middle = gene.start + (gene.end - gene.start[i])/2)) %>%
    left_join(markers %>%
                dplyr::select(marker.id, pos)) %>%
    mutate(gene.chr = chr,
            cis = TRUE)
  
  local_lod_table
}

######################################################################
###### 
###### Run local QTL (transcript to protein) through bmediatR
######
######################################################################
run_local_pqtl_bmediatR <- function (local_pqtl_table,
                                     protein_mat,
                                     transcript_mat,
                                     genoprobs,
                                     Z = NULL, Z_y = NULL, Z_M = NULL,
                                     w = NULL, w_y = NULL, w_M = NULL,
                                     phi_sq = c(1, 1, 1),
                                     use_emp_phi_sq = FALSE,
                                     ...) {
  
  partial_results <- complete_results <- reactive_results <- matrix(NA, nrow = nrow(local_pqtl_table), ncol = 12)
  rownames(partial_results) <- rownames(complete_results) <- rownames(reactive_results) <- local_pqtl_table$protein.id
  for (i in 1:nrow(local_pqtl_table)) {
    this_marker <- local_pqtl_table$marker.id[i]
    this_genoprobs <- qtl2::pull_genoprobpos(genoprobs, 
                                             marker = this_marker)
    # Empirical effect sizes
    if (use_emp_phi_sq) {
      emp_effect_size <- calc_trio_effect_sizes(y = protein_mat[,local_pqtl_table$protein.id[i]],
                                                m = transcript_mat[,local_pqtl_table$gene.id[i]],
                                                X = this_genoprobs,
                                                Z = Z, Z_y = Z_y, Z_m = Z_M,
                                                w = w, w_y = w_y, w_m = w_M,
                                                ...)
      phi_sq <- renormalize_effect_size_ratio(emp_effect_size)
    }
    # Partial
    partial <- bmediatR(y = protein_mat[,local_pqtl_table$protein.id[i]],
                        M = transcript_mat[,local_pqtl_table$gene.id[i]], 
                        X = this_genoprobs,
                        Z = Z, Z_y = Z_y, Z_M = Z_M,
                        w = w, w_y = w_y, w_M = w_M,
                        ln_prior_c = "partial",
                        phi_sq = phi_sq,
                        verbose = FALSE)
    partial_results[i,] <- partial$ln_post_c[1,]
    # Complete
    complete <- bmediatR(y = protein_mat[,local_pqtl_table$protein.id[i]],
                         M = transcript_mat[,local_pqtl_table$gene.id[i]], 
                         X = this_genoprobs,
                         Z = Z, Z_y = Z_y, Z_M = Z_M,
                         w = w, w_y = w_y, w_M = w_M,
                         ln_prior_c = "complete", 
                         phi_sq = phi_sq,
                         verbose = FALSE)
    complete_results[i,] <- complete$ln_post_c[1,]
    # Reactive
    reactive <- bmediatR(y = protein_mat[,local_pqtl_table$protein.id[i]],
                         M = transcript_mat[,local_pqtl_table$gene.id[i]], 
                         X = this_genoprobs,
                         Z = Z, Z_y = Z_y, Z_M = Z_M,
                         w = w, w_y = w_y, w_M = w_M,
                         ln_prior_c = "reactive", 
                         phi_sq = phi_sq,
                         verbose = FALSE)
    reactive_results[i,] <- reactive$ln_post_c[1,]
  }
  colnames(partial_results) <- colnames(complete_results) <- colnames(reactive_results) <- colnames(partial$ln_post_c)
  final_results <- list(partial = partial_results,
                        complete = complete_results,
                        reactive = reactive_results)
  final_results
}

######################################################################
###### 
###### Run local QTL (transcript to protein) LOD drop through
###### intermediate2
######
######################################################################
# Run local QTL through intermediate
run_local_pqtl_intermediate <- function (local_pqtl_table,
                                         protein_mat,
                                         transcript_mat,
                                         genoprobs,
                                         covar,
                                         ...) {
  
  lod_drop_results <- matrix(NA, nrow = nrow(local_pqtl_table), ncol = 4)
  rownames(lod_drop_results) <- local_pqtl_table$protein.id
  for (i in 1:nrow(local_pqtl_table)) {
    this_marker <- local_pqtl_table$marker.id[i]
    this_genoprobs <- qtl2::pull_genoprobpos(genoprobs, 
                                             marker = this_marker)
    this_M <- cbind(rep(1, nrow(transcript_mat)), transcript_mat[,local_pqtl_table$gene.id[i]])
    colnames(this_M) <- c("holder", "M")
    
    # LOD drop
    lod_drop <- intermediate2::mediation_scan(target = protein_mat[,local_pqtl_table$protein.id[i]],
                                              mediator = this_M,
                                              driver = this_genoprobs[rownames(protein_mat),], 
                                              covar = covar,
                                              annotation = data.frame(id = colnames(this_M),
                                                                      chr = c(1, 1), 
                                                                      pos = c(1, 2)))
    lod_drop_results[i,] <- c(lod_drop$lod[1], 
                              lod_drop$lod[2], 
                              lod_drop$lod[1] - lod_drop$lod[2],
                              (lod_drop$lod[1] - lod_drop$lod[2])/lod_drop$lod[1])
  }
  colnames(lod_drop_results) <- c("Y_lod", "Y_M_lod", "lod_drop", "lod_drop_prop")
  
  lod_drop_results
}


######################################################################
###### 
###### Plotting local QTL LOD score comparison
######
######################################################################
plot_local_pqtl_lod_comparison <- function (local_pqtl_table,
                                            bmediatR_results,
                                            model,
                                            col = c("seagreen4", "seagreen1", "skyblue", "goldenrod4", "goldenrod1", "gray"),
                                            include_count = TRUE,
                                            ncol = 3, nrow = 2) {
  
  ## Make data
  local_pqtl_data <- local_pqtl_table %>%
    left_join(exp(bmediatR_results[[model]]) %>%
                as.data.frame %>%
                rownames_to_column("protein.id"))
  # Grabbing most likely model
  local_pqtl_data$most_likely <- local_pqtl_data %>%
    dplyr::select(protein.id, contains(",")) %>%
    column_to_rownames("protein.id") %>%
    apply(1, function(x) names(x)[which.max(x)])
  # Grabbing maximum probability
  local_pqtl_data <- local_pqtl_data %>%
    mutate(max_prob = case_when(most_likely == "1,1,0" ~ `1,1,0`,
                                most_likely == "1,1,1" ~ `1,1,1`,
                                most_likely == "1,0,1" ~ `1,0,1`,
                                most_likely == "1,*,1" ~ `1,*,1`,
                                most_likely == "0,*,1" ~ `0,*,1`)) %>%
    rowwise %>%
    mutate(max_prob = ifelse(is.na(max_prob), 
                             sum(`0,0,0`, `0,1,0`, `1,0,0`, `0,0,1`, 
                                 `0,1,1`, `0,*,0`, `1,*,0`),
                             max_prob)) %>%
    as.data.frame
  
  # Simplify categories
  model_flag <- is.finite(bmediatR_results[[model]][1,])
  names(model_flag) <- colnames(bmediatR_results[[model]])
  ## Long names
  long_names <- c("other non-med", "other non-med", "other non-med", 
                  "complete med", "other non-med", "other non-med",
                  "co-local", "partial med", "other non-med",
                  "complete med (react)", "other non-med", "partial med (react)")
  names(long_names) <- colnames(bmediatR_results[[model]])
  
  col <- col[c(model_flag[c("1,1,0", "1,1,1", "1,0,1", "1,*,1", "0,*,1")], TRUE)]
  
  local_pqtl_data <- local_pqtl_data %>%
    mutate(best_model = case_when(most_likely == "1,1,1" ~ "partial med",
                                  most_likely == "1,1,0" ~ "complete med",
                                  most_likely == "1,0,1" ~ "co-local",
                                  most_likely == "1,*,1" ~ "partial med (react)",
                                  most_likely == "0,*,1" ~ "complete med (react)")) %>%
    mutate(best_model = ifelse(is.na(best_model), "other non-med", best_model)) %>%
    dplyr::mutate(best_model = factor(best_model, levels = c("complete med", "partial med", "co-local", "complete med (react)", "partial med (react)", "other non-med")))
  
  count_dat <- local_pqtl_data %>%
    dplyr::select(best_model) %>%
    group_by(best_model) %>%
    tally() %>%
    ungroup %>%
    mutate(pos_y = max(local_pqtl_data$protein_lod) * 0.9,
           pos_x = max(local_pqtl_data$transcript_lod) * 0.9)
  
  # Comparing lods by most likely model
  lods_p <- ggplot(data = local_pqtl_data, aes(x = transcript_lod, y = protein_lod, fill = best_model, size = max_prob)) + 
    geom_point(shape = 21, alpha = 0.7) + 
    geom_abline(intercept = 0, slope = 1, col = "black", linetype = "dashed") +
    coord_fixed()
  if (include_count) {
    lods_p <- lods_p +
      geom_text(data = count_dat, aes(label = n, x = pos_x, y = pos_y), col = "black", size = 4)
  }
  lods_p <- lods_p +
    scale_fill_manual(values = col) +
    scale_size(range = c(0.1, 3)) +
    xlab("eQTL LOD score") + ylab("pQTL LOD score") + 
    theme_bw() +
    theme(strip.background = element_rect(fill = "white"),
          strip.text.x = element_text(size = 12),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          legend.position = "bottom") +
    labs(size = "Highest\nposterior prob") +
    guides(fill = FALSE) +
    facet_wrap(~best_model, nrow = nrow, ncol = ncol)
  
  lods_p
}

######################################################################
###### 
###### Plotting local QTL effect correlation histograms
######
######################################################################
plot_local_pqtl_cor_comparison <- function (local_pqtl_table,
                                            bmediatR_results,
                                            model,
                                            col = c("seagreen4", "seagreen1", "skyblue", "goldenrod4", "goldenrod1", "gray"),
                                            ncol = 3, nrow = 2) {
  
  ## Make data
  # Calculate effect correlation
  cor_dat <- local_pqtl_table %>%
    dplyr::select(protein.id, gene.id, contains("protein"), contains("transcript"), matches("dist"))
  
  cor_dat$cor <- sapply(1:nrow(cor_dat), function(i) cor(as.numeric(cor_dat[i, paste("protein", LETTERS[1:8], sep = "_")]),
                                                         as.numeric(cor_dat[i, paste("transcript", LETTERS[1:8], sep = "_")])))
  
  ## Merging
  local_pqtl_data <- cor_dat %>%
    left_join(exp(bmediatR_results[[model]]) %>%
                as.data.frame %>%
                rownames_to_column("protein.id"))
  
  # Grabbing most likely model
  local_pqtl_data$most_likely <- local_pqtl_data %>%
    dplyr::select(protein.id, contains(",")) %>%
    column_to_rownames("protein.id") %>%
    apply(1, function(x) names(x)[which.max(x)])
  # Grabbing maximum probability
  local_pqtl_data <- local_pqtl_data %>%
    mutate(max_prob = case_when(most_likely == "1,1,0" ~ `1,1,0`,
                                most_likely == "1,1,1" ~ `1,1,1`,
                                most_likely == "1,0,1" ~ `1,0,1`,
                                most_likely == "1,*,1" ~ `1,*,1`,
                                most_likely == "0,*,1" ~ `0,*,1`)) %>%
    rowwise %>%
    mutate(max_prob = ifelse(is.na(max_prob), 
                             sum(`0,0,0`, `0,1,0`, `1,0,0`, `0,0,1`, 
                                 `0,1,1`, `0,*,0`, `1,*,0`),
                             max_prob)) %>%
    as.data.frame
  
  # Simplify categories
  model_flag <- is.finite(bmediatR_results[[model]][1,])
  names(model_flag) <- colnames(bmediatR_results[[model]])
  ## Long names
  long_names <- c("other non-med", "other non-med", "other non-med", 
                  "complete med", "other non-med", "other non-med",
                  "co-local", "partial med", "other non-med",
                  "complete med (react)", "other non-med", "partial med (react)")
  names(long_names) <- colnames(bmediatR_results[[model]])
  
  col <- col[c(model_flag[c("1,1,0", "1,1,1", "1,0,1", "1,*,1", "0,*,1")], TRUE)]
  
  local_pqtl_data <- local_pqtl_data %>%
    mutate(best_model = case_when(most_likely == "1,1,1" ~ "partial med",
                                  most_likely == "1,1,0" ~ "complete med",
                                  most_likely == "1,0,1" ~ "co-local",
                                  most_likely == "1,*,1" ~ "partial med (react)",
                                  most_likely == "0,*,1" ~ "complete med (react)")) %>%
    mutate(best_model = ifelse(is.na(best_model), "other non-med", best_model)) %>%
    dplyr::mutate(best_model = factor(best_model, levels = c("complete med", "partial med", "co-local", "complete med (react)", "partial med (react)", "other non-med")))
  
  # Comparing effect correlation by most likely model
  cor_p <- ggplot(data = local_pqtl_data, aes(x = cor, fill = best_model)) + 
    geom_histogram() + 
    scale_fill_manual(values = col) +
    xlab("Haplotype effect correlation") +
    theme_bw() +
    theme(strip.background = element_rect(fill = "white"),
          strip.text.x = element_text(size = 12),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10)) +
    guides(fill = FALSE) +
    facet_wrap(~best_model, nrow = nrow, ncol = ncol)
  
  cor_p
}

######################################################################
###### 
###### Plotting local QTL LOD drop histograms
######
######################################################################
plot_local_pqtl_lod_drop_comparison <- function (local_pqtl_table,
                                                 bmediatR_results,
                                                 lod_drop_results,
                                                 model,
                                                 col = c("seagreen4", "seagreen1", "skyblue", "goldenrod4", "goldenrod1", "gray"),
                                                 ncol = 3, nrow = 2) {
  
  ## Make data
  ## Merging
  local_pqtl_data <- local_pqtl_table %>%
    left_join(exp(bmediatR_results[[model]]) %>%
                as.data.frame %>%
                rownames_to_column("protein.id"))
  
  ## Merging in lod drop if provided
  local_pqtl_data <- local_pqtl_data %>%
    left_join(lod_drop_results %>%
                as.data.frame %>%
                rownames_to_column("protein.id"))
  
  # Grabbing most likely model
  local_pqtl_data$most_likely <- local_pqtl_data %>%
    dplyr::select(protein.id, contains(",")) %>%
    column_to_rownames("protein.id") %>%
    apply(1, function(x) names(x)[which.max(x)])
  # Grabbing maximum probability
  local_pqtl_data <- local_pqtl_data %>%
    mutate(max_prob = case_when(most_likely == "1,1,0" ~ `1,1,0`,
                                most_likely == "1,1,1" ~ `1,1,1`,
                                most_likely == "1,0,1" ~ `1,0,1`,
                                most_likely == "1,*,1" ~ `1,*,1`,
                                most_likely == "0,*,1" ~ `0,*,1`)) %>%
    rowwise %>%
    mutate(max_prob = ifelse(is.na(max_prob), 
                             sum(`0,0,0`, `0,1,0`, `1,0,0`, `0,0,1`, 
                                 `0,1,1`, `0,*,0`, `1,*,0`),
                             max_prob)) %>%
    as.data.frame
  
  # Simplify categories
  model_flag <- is.finite(bmediatR_results[[model]][1,])
  names(model_flag) <- colnames(bmediatR_results[[model]])
  ## Long names
  long_names <- c("other non-med", "other non-med", "other non-med", 
                  "complete med", "other non-med", "other non-med",
                  "co-local", "partial med", "other non-med",
                  "complete med (react)", "other non-med", "partial med (react)")
  names(long_names) <- colnames(bmediatR_results[[model]])
  
  col <- col[c(model_flag[c("1,1,0", "1,1,1", "1,0,1", "1,*,1", "0,*,1")], TRUE)]
  
  local_pqtl_data <- local_pqtl_data %>%
    mutate(best_model = case_when(most_likely == "1,1,1" ~ "partial med",
                                  most_likely == "1,1,0" ~ "complete med",
                                  most_likely == "1,0,1" ~ "co-local",
                                  most_likely == "1,*,1" ~ "partial med (react)",
                                  most_likely == "0,*,1" ~ "complete med (react)")) %>%
    mutate(best_model = ifelse(is.na(best_model), "other non-med", best_model)) %>%
    dplyr::mutate(best_model = factor(best_model, levels = c("complete med", "partial med", "co-local", "complete med (react)", "partial med (react)", "other non-med")))
  
  lod_drop_p <- ggplot(data = local_pqtl_data, aes(x = lod_drop_prop, fill = best_model)) + 
    geom_histogram() + 
    geom_vline(xintercept = 0, col = "black", linetype = "dashed") +
    scale_fill_manual(values = col) +
    xlab("LOD drop proportion") +
    theme_bw() +
    theme(strip.background = element_rect(fill = "white"),
          strip.text.x = element_text(size = 12),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10)) +
    guides(fill = FALSE) +
    facet_wrap(~best_model, nrow = nrow, ncol = ncol)
  
  lod_drop_p
}

######################################################################
###### 
###### Plotting allelic series spectrum
######
######################################################################
make_gradient <- function(left_spectrum = make_spectrum(c("cyan", "white"), n = 50),
                          right_spectrum = make_spectrum(c("white", "magenta"), n = 51),
                          mode = c("horizontal", "vertical"), 
                          title.cex = 1.5, 
                          axis.cex = 1.5, 
                          text.cex = 1.5,
                          line = 3) {
  mode <- mode[1]
  spectrum <- c(left_spectrum, right_spectrum[-1])
  if (mode == "horizontal") {
    barplot(height = rep(0.1, length(spectrum)), 
            width = 1/length(spectrum), 
            density = 1000, 
            angle = 90, 
            col = spectrum, 
            border = FALSE, 
            space = FALSE, 
            axes = FALSE)
    axis(1, at = c(0, 1), labels = c("low", "high"), cex.axis = axis.cex)
    mtext(side = 1, text = "Effect", line = line, 
          cex = text.cex)
  }
  else {
    barplot(height = rep(0.1, length(spectrum)), 
            width = 1/length(spectrum), 
            density = 1000, 
            angle = 0, 
            col = spectrum, 
            border = FALSE, 
            space = FALSE, 
            axes = FALSE, 
            horiz = TRUE)
    axis(4, las = 1, at = c(0, 1), labels = c("low", "high"), cex.axis = axis.cex)
    mtext(side = 1, text = "Effect", cex = title.cex, 
          line = line)
  }
}
make_spectrum <- function (col.range, 
                           col.spectrum, n = 1000) {
  
  if (is.null(col.range)) {
    if (col.spectrum == "gray") {
      spectrum <- gray(level = n:1/n)
    }
    else {
      spectrum <- do.call(what = eval(parse(text = paste0("colorRamps::", 
                                                          col.spectrum))), args = list(n = n))
    }
  }
  else {
    spectrum <- colorRampPalette(col.range)
    spectrum <- spectrum(n)
  }
  return(spectrum)
}

######################################################################
###### 
###### Function to process bmediatR grid plots (single row plots)
######
######################################################################
process_grid_plots <- function(ggplot,
                               strip_col = c("seagreen4", "seagreen1", "skyblue", "gray"),
                               text_col = c("white", "black", "black", "black"),
                               border_col = c(NA, NA, "black", NA),
                               border_lwd = c(0, 0, 4, 0)) {
  
  g <- ggplot_gtable(ggplot_build(ggplot))
  stripr <- which(grepl("strip-t", g$layout$name))
  k <- 1
  for (i in stripr) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- strip_col[k]
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$col <- border_col[k]
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$lwd <- border_lwd[k]
    
    j <- which(grepl("text", g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$children[[1]]$gp$col <- text_col[k]
    k <- k + 1
  }
  grid::grid.draw(g)
}

######################################################################
###### 
###### Function to process bmediatR grid plots (multiple row plots)
######
######################################################################
process_multi_row_grid_plots <- function(ggplot,
                                         strip_col = c("seagreen4", "seagreen1", "skyblue", "goldenrod4", "goldenrod1", "gray"),
                                         text_col = c("white", "black", "black", "white", "black", "black"),
                                         border_col = c(NA, NA, "black", NA, NA, NA),
                                         border_lwd = c(0, 0, 4, 0, NA, NA, NA)) {
  
  g <- ggplot_gtable(ggplot_build(ggplot))
  
  stripr <- which(grepl("strip-t", g$layout$name))
  k <- 1
  for (i in stripr) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- strip_col[k]
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$col <- border_col[k]
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$lwd <- border_lwd[k]
    
    j <- which(grepl("text", g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$children[[1]]$gp$col <- text_col[k]
    k <- k + 1
  }
  grid::grid.draw(g)
}

######################################################################
###### 
###### Plot posterior probabilities as boxplots
######
######################################################################
multi_sim_boxplot <- function(post_odds_mat,
                              include_lines = TRUE,
                              box_col = c("seagreen4", "seagreen1", "goldenrod1", "goldenrod4", "skyblue", "gray")) {
  
  model_flag <- !as.character(post_odds_mat[1,]) == "-Inf"
  names(model_flag) <- colnames(post_odds_mat)
  
  ## Long names
  long_names <- c("other non-med", "other non-med", "other non-med", 
                  "complete med", "other non-med", "other non-med",
                  "co-local", "partial med", "other non-med",
                  "complete med (react)", "other non-med", "partial med (react)")
  names(long_names) <- names(model_flag)
  
  box_col <- box_col[c(model_flag[c("1,1,0", "1,1,1", "1,*,1", "0,*,1", "1,0,1")], TRUE)]
  
  post_odds_dat <- post_odds_mat %>%
    as.data.frame %>%
    rowid_to_column %>%
    mutate("partial med" = exp(`1,1,1`),
           "complete med" = exp(`1,1,0`),
           "co-local" = exp(`1,0,1`),
           "partial med (react)" = exp(`1,*,1`),
           "complete med (react)" = exp(`0,*,1`)) %>%
    dplyr::select(-c("1,1,1", "1,1,0", "1,0,1", "1,*,1", "0,*,1"))
  post_odds_dat <- post_odds_dat %>%
    dplyr::select(-contains(",")) %>%
    bind_cols(post_odds_dat %>%
                dplyr::select(contains(",")) %>%
                transmute("other non-med" = rowSums(exp(.)))) %>%
    gather(key = "model", value = "prob", -rowid)
  
  ## Set factors
  models_use <- unique(long_names[model_flag])
  post_odds_dat <- post_odds_dat %>%
    filter(model %in% models_use) %>%
    mutate(model = factor(model, levels = c("complete med", "partial med", "partial med (react)", 
                                            "complete med (react)", "co-local", "other non-med")))
  
  box_theme <- theme(axis.line = element_line(colour = "black"),
                     plot.title = element_text(hjust = 0.5, size = 16, face = "plain"), 
                     axis.title.x = element_blank(),
                     axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "plain"),
                     axis.title.y = element_text(size = 14, face = "plain"),
                     axis.text.y = element_text(size = 14, face = "plain"),
                     axis.ticks.y = element_blank(),
                     legend.title = element_text(size = 14),
                     legend.text = element_text(size = 14))
  
  if (include_lines) {
    ggplot(data = post_odds_dat, aes(y = prob, x = model, color = model)) + 
      geom_line(aes(group = rowid), color = "gray80", alpha = 0.1) + 
      geom_boxplot() +
      scale_color_manual(values = box_col) +
      theme_bw() + ylim(0, 1) +
      ylab("Posterior model probability") + xlab("Model") + box_theme + guides(color = "none")
  } else {
    ggplot(data = post_odds_dat, aes(y = prob, x = model, color = model)) + 
      geom_boxplot() +
      scale_color_manual(values = box_col) +
      theme_bw() + ylim(0, 1) +
      ylab("Posterior model probability") + xlab("Model") + box_theme + guides(color = "none")
  }
}

