## Simple muller plot with drug scheduling annotations

plot_muller <- function(df_pop, sample, rep) {
  
  
  library(ggmuller)
  library(tidyverse)
  library(zoo)
  library(ggthemes)
  
  MFREQ <- 0.001
  
  label_at <- function(n) function(x) ifelse(x %% n == 0, x, "")
  
  line_2d <- case_when(
    # grepl("Trametinib ->", sample) ~ 74,
    # grepl("->Trametinib", sample) ~ 74,
    TRUE ~ 5)
  
  no_breaks <- case_when(
    grepl("Trametinib", sample) ~ 3,
    TRUE ~ 2)
  
  line_exp <- case_when(
    
    grepl("Trametinib -> ",sample) ~ 126, 
    grepl("-> Trametinib", sample) ~ 126,
    TRUE ~ 47)
  
  savedir <- "../res/muller_plots_erc/"
  
  df_pop <- df_pop %>% group_by(Identity) %>% mutate(m_freq = max(Frequency))
  sample_init <- "nnn"
  if(grepl("->",sample )){
    sample_init <- strsplit(sample, split = "->")[[1]][1] %>% trimws()
    
    line_exp_pre <- case_when(
      
      grepl("Trametinib -> ",sample_init) ~ 126, 
      grepl("-> Trametinib", sample_init) ~ 126,
      TRUE ~ 47)
  } 
  
  line_nd <- case_when(
    
    grepl("Trametinib ->",sample) ~ 87, 
    grepl("-> Trametinib", sample) ~ 69,
    TRUE ~ 3)
  
  data_pop <- df_pop %>% filter((sample == !!sample) | (sample == !!sample_init)) %>% filter( rep == !!rep, m_freq > !!MFREQ)  %>% select(-sample, -rep, -Frequency, -m_freq) 
  data_edges <- data.frame(Parent = "0", Identity = data_pop$Identity) %>% unique()
  data_pop <- rbind(data_pop[1,,drop = F] %>% mutate(Identity = "0", Population = 100, Time = 0), data_pop)
  
  
  Muller_df <- get_Muller_df(data_edges, data_pop, cutoff = 0, start_positions = 0.5)
  
  
  palette <- readRDS("~/data/sc_project/javi_sc_rna_seq/intermediate_results//palette_enriched_100.rds")
  no_in_palette <- setdiff(Muller_df$Identity %>% unique(), names(palette))
  new_pal <- rep("grey50", times = length(no_in_palette))
  names(new_pal) <- no_in_palette
  palette <- c(new_pal, palette)
  
  p_mull_smoo <- Muller_plot(Muller_df %>% mutate(Time = Time - line_2d)  %>% group_by(Identity) %>% 
                               mutate(Frequency = zoo::rollmean(Frequency, 11, 
                                                                fill = c(Frequency[1],NA,
                                                                         Frequency[length(Frequency)]),
                                                                align = "center"))
                             , conceal_edges = F) + scale_fill_manual(values = palette, na.value = "white") +
    coord_cartesian(clip="off", xlim = c(line_2d, max(data_pop$Time))) + theme_tufte()+ theme(legend.position = "none", text=element_text(size=FONT_SIZE,  family="Arial"),
                                                                                              axis.text.x =element_text(size=8.5,  family="Arial"),
                                                                                              axis.text.y =element_text(size=8.5,  family="Arial"),
                                                                                              axis.ticks = element_line(size = 0.1),
                                                                                              plot.margin = margin(0,0,0,0, 'cm')) +
    scale_x_continuous(breaks = unique(data_pop$Time)[-c(1:no_breaks)] - 7, limits = c(0, max(data_pop$Time))) +
    
    xlab("Time (days)") + 
    scale_y_continuous(labels = 25 * (0:4), name = "Barcode frequency") + Seurat::RotatedAxis()
  
  # if(sample_init != "nnn") {
  #   p_mull_smoo <- p_mull_smoo + geom_vline(xintercept = as.numeric(line_nd), color = "darkred", linetype = "dashed") + 
  #     geom_vline(xintercept = as.numeric(line_exp), color = "black", linetype = "dashed") + geom_vline(xintercept = as.numeric(line_exp_pre), color = "black", linetype = "dashed")
  # } else {
  #   p_mull_smoo <- p_mull_smoo + geom_vline(xintercept = as.numeric(line_exp), color = "black", linetype = "dashed")
  # }
  
  drug1 <- strsplit(sample, split = " -> ")[[1]][1]
  drug2 <- strsplit(sample, split = " -> ")[[1]][2]
  
  df_annot <- data.frame(x = c(0, line_exp_pre, line_nd, line_exp), stage = c(paste0(drug1, " treatment"), "Re-growth", paste0(drug2, " treatment"), "Re-growth"), w = c(
    line_exp_pre - 0, line_nd - line_exp_pre, line_exp - line_nd,  p_mull_smoo$data$Time %>% max() - line_exp)
  )
  
  cols_annot <- c(drug_cols[drug1], "grey80", drug_cols[drug2], "grey80")
  names(cols_annot) <- c(paste0(drug1, " treatment"), "Re-growth", paste0(drug2, " treatment"), "Re-growth")
  
  annot_barplot <- ggplot(df_annot, aes(xmin = x, xmax = x + w, fill = stage, ymin = 0, ymax = 1)) +
    ggtitle(paste0(sample, " replicate ", rep)) + 
    geom_rect(color = "grey10")+ theme_void() + scale_fill_manual("",values = cols_annot)+
    theme(legend.position = "top", text=element_text(size=FONT_SIZE,  family="Arial"), plot.margin = margin(0,0,0,0, 'cm')) + coord_cartesian(clip="off", xlim = c(line_2d, max(data_pop$Time)))+
    scale_x_continuous(breaks = unique(data_pop$Time)[-c(1:no_breaks)] - 7, limits = c(0, max(data_pop$Time)))
  
  p_final <- cowplot::plot_grid(annot_barplot, p_mull_smoo, align = "v",axis = "l", ncol = 1, rel_heights = c(0.3,1))
  return(p_final)
  
}



plot_muller2 <- function(df_pop, sample, rep, palette, plot_bar = T, window_size = 9, index_to_rm = 3) {
  
  
  library(ggmuller)
  library(tidyverse)
  library(zoo)
  library(ggthemes)
  
  MFREQ <- 0.001
  
  label_at <- function(n) function(x) ifelse(x %% n == 0, x, "")
  
  
  if(grepl("Oxa", sample)){
    line_exp <- OXA_exp[[sample]]
  } else {
    line_exp <- ERK_exp[[sample]]
  }
  
  
  df_pop <- df_pop %>% dplyr::group_by(Identity) %>% dplyr::mutate(m_freq = max(Frequency)) %>% 
    ungroup()
  sample_init <- "nnn"
  
  
  data_pop <- df_pop %>% dplyr::filter((sample == !!sample) | (sample == !!sample_init)) %>% dplyr::filter( rep == !!rep, m_freq > !!MFREQ)  %>% 
    dplyr::select(-sample, -rep, -Frequency, -m_freq) 
  data_edges <- data.frame(Parent = "0", Identity = data_pop$Identity %>% unique()) 
  data_pop <- rbind(data_pop[1,,drop = F] %>% 
                      dplyr::mutate(Identity = "0", Population = 100, Time = 0), data_pop)
  
  Muller_df <- get_Muller_df(data_edges, data_pop, cutoff = 0, start_positions = 0.5, smooth_start_points = TRUE)
  

  no_in_palette <- setdiff(Muller_df$Identity %>% unique(), names(palette))
  new_pal <- rep("grey50", times = length(no_in_palette))
  names(new_pal) <- no_in_palette
  palette <- c(new_pal, palette)
  
  p_mull_smoo <- Muller_plot(Muller_df %>% dplyr::mutate(Time = as.double(Time))  %>% dplyr::group_by(Group_id) %>% 
                               arrange(Group_id, Time, row_number()) %>%
                               dplyr::mutate(Frequency = zoo::rollmean(Frequency, window_size, 
                                                                       fill = c(Frequency[1],NA,
                                                                                Frequency[length(Frequency)]), align = "center")) %>%  ungroup()
                             , conceal_edges = F) + scale_fill_manual(values = palette, na.value = "white") +
    coord_cartesian(clip="off", xlim = c(gtools::mixedsort(unique(data_pop$Time))[index_to_rm], max(data_pop$Time))) + theme_tufte()+ theme(legend.position = "none", text=element_text(size=FONT_SIZE,  family="Arial"),
                                                                                                                                            axis.text.x =element_text(size=8.5,  family="Arial"),
                                                                                                                                            axis.text.y =element_text(size=8.5,  family="Arial"),
                                                                                                                                            axis.ticks = element_line(size = 0.1),
                                                                                                                                            plot.margin = margin(0,0,0,0, 'cm')) +
    scale_x_continuous(breaks = gtools::mixedsort(unique(data_pop$Time))[-c(1:(index_to_rm - 1))], labels = (gtools::mixedsort(unique(data_pop$Time)) - gtools::mixedsort(unique(data_pop$Time))[index_to_rm])[-c(1:(index_to_rm - 1))],
                       limits = c(gtools::mixedsort(unique(data_pop$Time))[index_to_rm], max(data_pop$Time))) +
    
    xlab("Time (days)") + 
    #scale_y_continuous(labels = 25 * (0:4), name = "Barcode frequency") +
    Seurat::RotatedAxis()
  
  # if(sample_init != "nnn") {
  #   p_mull_smoo <- p_mull_smoo + geom_vline(xintercept = as.numeric(line_nd), color = "darkred", linetype = "dashed") + 
  #     geom_vline(xintercept = as.numeric(line_exp), color = "black", linetype = "dashed") + geom_vline(xintercept = as.numeric(line_exp_pre), color = "black", linetype = "dashed")
  # } else {
  #   p_mull_smoo <- p_mull_smoo + geom_vline(xintercept = as.numeric(line_exp), color = "black", linetype = "dashed")
  # }
  if(plot_bar) {
    
    drug <- strsplit(sample, split = "-")[[1]][3]  
    
    df_annot <- data.frame(x = c(0, line_exp), stage = c(paste0(drug, " treatment"), "Re-growth"), w = c(
      line_exp - 0,  max(data_pop$Time) - line_exp))
    df_annot$x <- df_annot$x + gtools::mixedsort(unique(data_pop$Time))[index_to_rm]
    df_annot$w <- ifelse(df_annot$x + df_annot$w > max(unique(data_pop$Time)),max(unique(data_pop$Time)) - df_annot$x,df_annot$w  )
    
    cols_annot <- c(drug_cols[drug], "grey80")
    names(cols_annot) <- c(paste0(drug, " treatment"), "Re-growth")
    
    annot_barplot <- ggplot(df_annot, aes(xmin = x, xmax = x + w, fill = stage, ymin = 0, ymax = 1)) +
      ggtitle(paste0(sample, " replicate ", rep)) + 
      geom_rect(color = "grey10")+ theme_void() + scale_fill_manual("",values = cols_annot)+
      theme(legend.position = "top", text=element_text(size=FONT_SIZE,  family="Arial"), plot.margin = margin(0,0,0,0, 'cm')) + coord_cartesian(clip="off", xlim = c(gtools::mixedsort(unique(data_pop$Time))[index_to_rm], max(data_pop$Time)))+
      scale_x_continuous(limits = c(gtools::mixedsort(unique(data_pop$Time))[index_to_rm], max(data_pop$Time)))
    
    p_final <- cowplot::plot_grid(annot_barplot, p_mull_smoo, align = "v",axis = "l", ncol = 1, rel_heights = c(0.3,1))
    return(p_final)
  } else {
    return(p_mull_smoo)
  }
  
  
}




from_bc_to_barplot_colors <- function(seur_obj, bc, name, clr){
  library(dplyr)
  meta <- seurat_scrna_seq_has_bc@meta.data
  name_str = glue::glue(" <span style = 'color:{clr};'>{name}</span> ")
  new_col <- meta %>% mutate( new_col = case_when(
    bc44 == !!bc & stage == "parental" ~ paste0("Parental (", name_str, ")"),
    bc44 == !!bc & stage == "under_drug" ~ paste0("Under Drug (", name_str, ")"),
    bc44 == !!bc & stage == "expansion" ~ paste0("Re-Growth (", name_str, ")"),
    bc44 != !!bc & stage == "parental" ~ "Parental (others)"
  )) %>% pull(new_col)
  new_col <- factor(new_col, levels = c("Parental (others)", paste0("Parental (", name_str, ")"),
                                        paste0("Under Drug (", name_str, ")"), 
                                        paste0("Re-Growth (", name_str, ")")))
}


from_bc_to_barplot_colors_batch_2 <- function(meta, bc, name, clr){
  library(dplyr)
  name_str = glue::glue(" <span style = 'color:{clr};'>{name}</span> ")
  new_col <- meta %>% mutate( new_col = case_when(
    bc44 == !!bc & stage == "DMSO" ~ paste0("DMSO (", name_str, ")"),
    bc44 == !!bc & stage == "under_drug" ~ paste0("Under Drug (", name_str, ")"),
    bc44 == !!bc & stage == "expansion" ~ paste0("Re-Growth (", name_str, ")"),
    bc44 != !!bc & stage == "DMSO" ~ "DMSO (others)"
  )) %>% pull(new_col)
  new_col <- factor(new_col, levels = c("DMSO (others)", paste0("DMSO (", name_str, ")"),
                                        paste0("Under Drug (", name_str, ")"), 
                                        paste0("Re-Growth (", name_str, ")")))
}


from_bc_to_barplot_colors_batch_2_archr <- function(meta, bc, name, clr){
  library(dplyr)
  name_str = glue::glue(" <span style = 'color:{clr};'>{name}</span> ")
  new_col <- meta %>% mutate( new_col = case_when(
    bc44 == !!bc & stage == "DMSO" ~ paste0("DMSO (", name_str, ")"),
    bc44 == !!bc & stage == "expansion" ~ paste0("Re-Growth (", name_str, ")"),
    bc44 != !!bc & stage == "DMSO" ~ "DMSO (others)"
  )) %>% pull(new_col)
  new_col <- factor(new_col, levels = c("DMSO (others)", paste0("DMSO (", name_str, ")"),
                                        paste0("Re-Growth (", name_str, ")")))
}

from_bc_to_barplot_colors_archr <- function(meta, bc, name, clr){
  library(dplyr)
  name_str = glue::glue(" <span style = 'color:{clr};'>{name}</span> ")
  new_col <- meta %>% mutate( new_col = case_when(
    bc44 == !!bc & stage == "parental" ~ paste0("Parental (", name_str, ")"),
    bc44 == !!bc & stage == "expansion" ~ paste0("Re-Growth (", name_str, ")"),
    bc44 != !!bc & stage == "parental" ~ "Parental (others)"
  )) %>% pull(new_col)
  new_col <- factor(new_col, levels = c("Parental (others)", paste0("Parental (", name_str, ")"),
                                        paste0("Re-Growth (", name_str, ")")))
}


plot_archetypes_simplex <- function(res, distance_type = "euclidean", color_by = NULL, 
                                    subsample = NULL, point_size = 1, 
                                    label_size = 12, label_title = "Group", 
                                    color_palette = NULL, 
                                    group_by_var = NULL, group_by_var_2 = NULL,
                                    filter_cells = NULL, legend = FALSE, title_to_add = "") {
  
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(geosphere)
  library(ggtext)
  library(RColorBrewer)
  
  # Extract archetypes or subset if specified
  arcs <- if (!is.null(subsample)) res$A[subsample, ] else res$A
  
  if(!is.null(filter_cells)){
    arcs <- arcs[filter_cells,]
    color_by <- color_by[filter_cells]
    if(!is.null(group_by_var)) group_by_var <- group_by_var[filter_cells]
    if(!is.null(group_by_var_2)) group_by_var_2 <- group_by_var_2[filter_cells]
  }
  
  if(!is.null(color_by)){
    arcs <- arcs[!is.na(color_by),]
    if(!is.null(group_by_var)) group_by_var <- group_by_var[!is.na(color_by)]
    if(!is.null(group_by_var_2)) group_by_var_2 <- group_by_var_2[!is.na(color_by)]
    color_by <- color_by[!is.na(color_by)]
  }
  
  
  # Calculate distance matrix and solve TSP
  dist_matrix <- dist(t(arcs), method = distance_type, diag = T, upper = T)
  library(TSP)
  tsp_solution <- solve_TSP(TSP(dist_matrix))
  order <- tsp_solution %>% as.integer()
  dist_matrix <- dist_matrix %>% as.matrix()
  
  # Calculate angles and radii for the polar plot
  dists <- sapply(1:(length(order)-1), function(i) dist_matrix[order[i], order[i+1]])
  dists = c(dists, dist_matrix[order[length(order)],order[1]])
  cumulative_distances <- c(0, cumsum(dists / sum(dists)))
  angles <- 360 * cumulative_distances[-1]
  angles_rad <- pracma::deg2rad(angles)
  
  # Calculate Cartesian coordinates for the plot
  x <- rowSums(apply(arcs,1, function(x) x* cos(angles_rad)) %>% t)
  y <- rowSums(apply(arcs,1, function(x) x* sin(angles_rad)) %>% t)
  data <- data.frame(Theta = atan2(y, x), Radius = sqrt(x^2 + y^2))
  #data$Theta <- ifelse(data$Theta < 0, 2*pi - data$Theta, data$Theta)
  
  # If color_by is specified, use it
  if (is.null(color_by)) {
    data$Color <- data$Theta
    colfunc <- grDevices::colorRampPalette(brewer.pal(name = "Set1",ncol(arcs))) # the choice of colors is yours!
    data$Color <- colfunc(data$Theta)
  } else {
    data$Color <- color_by 
  }
  if(!is.null(group_by_var))
    data$gp = group_by_var
  
  if(!is.null(group_by_var_2))
    data$gp2 = group_by_var_2
  
  # Plot using ggplot2
  p <- ggplot(data, aes(x = Theta, y = Radius, color = Color, alpha = Radius)) +
    geom_point(size = point_size, position=position_jitter(h=0.01, w=0.01)) +
    coord_polar() +
    scale_x_continuous(limits = c(-pi,  pi), labels = paste0("arc", 1:ncol(arcs))[order], breaks = atan2(sin(angles_rad), cos(angles_rad))) +
    theme_minimal() +
    labs(title = paste0("Archetype projection in a 2D polytope", title_to_add), color = label_title) + guides(alpha="none") +
    xlab("") + ylab("")  + theme(axis.text.y = element_blank(), panel.grid.minor = element_blank()) + ylim(0,0.9)
  # Add labels if needed
  if (is.null(color_by)) {
    p <- p + theme(legend.position = "none") 
  } else {
    
    p <- p + theme(legend.title = element_text(size = label_size), text = element_text(family = "Arial"),
                   legend.text = element_markdown(family = "Arial"), strip.text = element_markdown(family = "Arial")) + 
      scale_color_manual(values = color_palette)
  }
  if(!is.null(group_by_var) & is.null(group_by_var_2))
    p <- p + facet_grid(~gp)
  if(!is.null(group_by_var) & !is.null(group_by_var_2))
    p <- p + facet_grid(gp2~gp)
  if(!legend){
    p <- p + theme(legend.position = "none") 
  } else{
    p <- p + theme(legend.position = "bottom") 
  }
  p
}



archetype_barplot <- function(geno_df, palette, col_to_group, min_n = 10) {
  palette_archs <- RColorBrewer::brewer.pal(length(geno_df$archs %>% unique), "Set3")
  df_plot <- geno_df 
  
  p1 <- ggplot(df_plot %>% group_by(!!as.name(col_to_group),stage , drug) %>% mutate(N = length(unique(key)), group = paste0(!!as.name(col_to_group), " (", N ,")") ) %>% ungroup() %>% 
                 select(-key, -weights, -!!as.name(col_to_group)) %>% unique() %>% filter(N > !!min_n) , aes(x =group , y = avg_prob, fill = archs )) + 
    geom_col() + theme_bw() + labs(xlab = "%") + facet_wrap(stage ~ drug, scales = "free") + RotatedAxis()   + scale_fill_manual("Archs", values  = palette_archs)
  
  p2 <- ggplot(df_plot %>% group_by(!!as.name(col_to_group), stage) %>% mutate(N = length(unique(key)), group = paste0(!!as.name(col_to_group), " (", N ,")") ) %>% ungroup() %>% 
                 select(-key, -weights, -!!as.name(col_to_group)) %>% unique()  %>% filter(N > !!min_n), aes(x = group, y = avg_prob, fill = archs )) + 
    geom_col() + theme_bw() + labs(xlab = "%") + facet_wrap(stage ~ ., scales = "free") + RotatedAxis()   + scale_fill_manual("Archs", values  = palette_archs)
  p3 <- ggplot(df_plot %>% group_by(!!as.name(col_to_group)) %>% mutate(N = length(unique(key)), group = paste0(!!as.name(col_to_group), " (", N ,")") ) %>% ungroup() %>% 
                 select(-key, -weights, -!!as.name(col_to_group), drug) %>% unique()  %>% filter(N > !!min_n) , aes(x =group, y = avg_prob, fill = archs )) + 
    geom_col() + theme_bw() + labs(xlab = "%") + facet_wrap(drug ~ ., scales = "free") + RotatedAxis()   + scale_fill_manual("Archs", values  = palette_archs)
  p4 <- ggplot(df_plot %>% group_by(!!as.name(col_to_group)) %>% mutate(N = length(unique(key)), group = paste0(!!as.name(col_to_group), " (", N ,")") ) %>% ungroup() %>% 
                 select(-key, -weights, -!!as.name(col_to_group)) %>% unique()  %>% filter(N > !!min_n), aes(x =group, y = avg_prob, fill = archs )) + 
    geom_col() + theme_bw() + labs(xlab = "%")  + RotatedAxis()   + scale_fill_manual("Archs", values  = palette_archs)
  # p2 <-  ggplot(geno_df %>% select(-key, -weights) %>% unique(), aes(x = archs, y = avg_prob, fill = !!as.name(col_to_group) )) + 
  #   geom_col() + theme_bw() + labs(xlab = "N") + facet_wrap(stage ~ drug, scales = "free") + RotatedAxis()   + scale_fill_manual("Clusters", values  = palette)
  
  return(list("by_all" = p1,"by_stage" = p2, "by_drug" = p3, "by_none" = p4))
}



archetype_waddington_1d <- function(arch_table, drug, stage, palette , group_by_col = "bc44",
                                    greater_than = NULL, return_bottom = TRUE, include_na = FALSE, return_top = FALSE, reverse_order = FALSE) {
  
  library(tidyverse)
  library(ggpubr)
  
  arch_table_filt <- arch_table %>% dplyr::filter(drug == !!drug, stage == !!stage) 
  
  if(!include_na) arch_table_filt <- arch_table_filt %>% dplyr::filter(!is.na(!!group_by_col)) 
  
  n_arch <- length(unique(arch_table_filt$archs))
  if(is.null(greater_than)) {
    all_barcodes <- unique(arch_table_filt %>% dplyr::pull(!!group_by_col))
  } else {
    tb_bc <- arch_table_filt %>% select(key, !!group_by_col) %>% unique() %>% dplyr::pull(!!group_by_col) %>% table()
    all_barcodes <- names(tb_bc)[tb_bc > greater_than]
  }
  
  plots <- vector("list", length(all_barcodes))
  names(plots) <- all_barcodes
  for(i in all_barcodes){
    plots[[i]] <- archetype_waddington_1d_aux(arch_table_filt, drug, stage, i, n_arch, palette, return_bottom, group_by_col, return_top, reverse_order)
  }
  
  return(plots)
}

archetype_waddington_1d_aux <- function(arch_table_filt, drug, stage, bc44_inp, n_arch, palette, return_bottom, group_by_col, return_top, reverse_order = FALSE) {
  
  
  arch_names <- arch_table_filt$archs %>% unique()
  
  x <- seq(0, 2 * pi, length.out = 10000)
  probs_arc <- arch_table_filt %>% dplyr::ungroup() %>% 
    dplyr::mutate(bc44_plot = dplyr::if_else(!!as.name(group_by_col) != !!bc44_inp | is.na(!!group_by_col),  "others", bc44_inp)) %>% 
    dplyr::group_by(bc44_plot, archs) %>% dplyr::summarize(avg_prob = mean(weights)) %>% 
    dplyr::ungroup() %>% dplyr::group_by(bc44_plot) %>% dplyr::mutate(avg_prob_sum = sum(avg_prob)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(avg_prob = avg_prob / avg_prob_sum) %>% 
    dplyr::select(archs, avg_prob, bc44_plot) %>% unique() %>% dplyr::arrange(archs) 
  
  n_barcode <- arch_table_filt %>% dplyr::filter(!!as.name(group_by_col) == !!bc44_inp) %>% pull(key) %>% unique() %>% length()
  n_others <- arch_table_filt %>% dplyr::filter(!!as.name(group_by_col) != !!bc44_inp | is.na(!!group_by_col) ) %>% dplyr::pull(key) %>% unique() %>% length()
  probs_arc_other <- probs_arc %>% dplyr::filter(bc44_plot == "others") %>% dplyr::pull(avg_prob)
  probs_arc <-  probs_arc %>% dplyr::filter(bc44_plot == !!bc44_inp) %>% pull(avg_prob)
  names(probs_arc) <- arch_names
  names(probs_arc_other) <- arch_names
  y <-  cos(x * n_arch) - 1
  y_others <- y
  lims = cumsum(2 * pi / n_arch * rep(1, n_arch))
  lims = c(0, lims)
  
  x_labels = 2 * pi / (n_arch * 4) * seq(0,n_arch * 4 - 1, length.out = n_arch * 4)
  x_labels <- x_labels[seq(3,n_arch * 4, by = 4)]
  
  y_labels <- cos(x_labels * n_arch) - 1
  y_others_labels <- y_labels
  
  for(i in seq(1, length(lims) - 1)) {
    y[x > lims[i] & x <= lims[i+1]] = y[x > lims[i] & x <= lims[i+1]] * probs_arc[i]
    y_labels[x_labels > lims[i] & x_labels <= lims[i+1]] =  y_labels[x_labels > lims[i] & x_labels <= lims[i+1]] * probs_arc[i]
    
    y_others[x > lims[i] & x <= lims[i+1]] = y_others[x > lims[i] & x <= lims[i+1]] * probs_arc_other[i]
    y_others_labels[x_labels > lims[i] & x_labels <= lims[i+1]] =  y_others_labels[x_labels > lims[i] & x_labels <= lims[i+1]] * probs_arc_other[i]
    
  }
  
  df_plot_density_bc <- data.frame(x,y, bc = bc44_inp ) 
  df_plot_density_others <- data.frame(x,y = y_others, bc = "others" ) 
  df_plot_density <- rbind(df_plot_density_bc, df_plot_density_others)
  
  df_plot_labels_bc <- data.frame(x = x_labels, y =  y_labels - 0.03, lab = round(probs_arc, 2), bc = bc44_inp, 
                                  arch = arch_names)
  df_plot_labels_other <- data.frame(x = x_labels, y =  y_others_labels - 0.03, 
                                     lab = round(probs_arc_other, 2), bc = "others", arch = arch_names)
  df_plot_labels <- rbind(df_plot_labels_bc, df_plot_labels_other)
  
  
  
  
  if(return_top) {
    p_density <- ggplot(df_plot_density %>% filter(bc != "others"), aes(x, y, color = bc)) + geom_line(linewidth = 1.5) +
      geom_label(df_plot_labels  %>% filter(bc != "others"), mapping = aes(label = lab, color = bc)) +
      geom_label(df_plot_labels  %>% filter(bc != "others"), mapping = aes(y = 0.05,label = arch), color = "grey20") +
      scale_color_manual(values = palette) +
      ggridges::theme_ridges()  + Seurat::RotatedAxis() + scale_x_continuous(breaks = x_labels,labels = arch_names)  +
      theme(axis.text.y = element_blank(), panel.grid = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),
            legend.position = "None", strip.text = element_blank(), panel.spacing = unit(0, "lines"))+ ylab("") + xlab("") +
      ggtitle(glue::glue("Archetypal landscape for {drug} in {stage} population"), 
              subtitle = glue::glue("Group {bc44_inp} ({n_barcode}) ")) 
  } else {
    p_density <- ggplot(df_plot_density, aes(x, y, color = bc)) + geom_line(linewidth = 1.5) +
      geom_label(df_plot_labels, mapping = aes(label = lab, color = bc)) +
      geom_label(df_plot_labels, mapping = aes(y = 0.05,label = arch), color = "grey20") +
      scale_color_manual(values = palette) +
      ggridges::theme_ridges()  + Seurat::RotatedAxis() + scale_x_continuous(breaks = x_labels,labels = arch_names)  +
      theme(axis.text.y = element_blank(), panel.grid = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),
            legend.position = "None", strip.text = element_blank(), panel.spacing = unit(0, "lines"))+ ylab("") + xlab("") +
      ggtitle(glue::glue("Archetypal landscape for {drug} in {stage} population"), 
              subtitle = glue::glue("Group {bc44_inp} ({n_barcode}) vs others ({n_others})")) 
    if(reverse_order){
      p_density <- p_density + facet_wrap(~ forcats::fct_rev(bc), ncol = 1)
    }  else {
      p_density <- p_density + facet_wrap(~ bc, ncol = 1)
    }
    
  }
  
  
  if(!return_bottom) return(p_density)
  
  # p_density <- ggplot(df_plot_density, aes(x, y, color = bc)) + geom_line(linewidth = 1.5) + 
  #   geom_label(df_plot_labels, mapping = aes(label = lab, color = bc)) + scale_color_manual(values = palette) +
  #   ggridges::theme_ridges()  + Seurat::RotatedAxis() + scale_x_continuous(breaks = x_labels,labels = arch_names)  +
  #   theme(axis.text.y = element_blank(), panel.grid = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),
  #         legend.position = "None", strip.text = element_blank(), panel.spacing = unit(0, "lines"))+ ylab("") + xlab("") +
  #   ggtitle("") + facet_wrap(~ bc, ncol = 1)
  
  arch_table_filt_boxplot <- arch_table_filt %>% ungroup() %>% 
    dplyr::mutate(bc44_plot = if_else(!!as.name(group_by_col) != !!bc44_inp | is.na(!!group_by_col),  "others", bc44_inp))
  
  palette_box <- c(palette[bc44_inp], "grey50")
  names(palette_box) <- c(bc44_inp, "others")
  
  p_boxplot <- ggplot(arch_table_filt_boxplot, aes(x = archs, y = weights, fill = bc44_plot)) + geom_boxplot() + 
    ggpubr::stat_compare_means(aes(label = ..p.signif..), method = "wilcox.test") + ggridges::theme_ridges() + xlab("") +
    scale_x_discrete(labels = arch_names)  + ylab("Archetype weight") + 
    scale_fill_manual("Group", values = palette_box) + 
    theme(panel.grid = element_blank(), legend.position = "bottom")  + 
    Seurat::RotatedAxis() 
  
  cowplot::plot_grid(p_density, p_boxplot, ncol = 1, align = "v", axis = "brtl", rel_heights = c(1.,1))
}




archetype_barplot <- function(geno_df, palette, col_to_group, min_n = 10) {
  palette_archs <- RColorBrewer::brewer.pal(length(geno_df$archs %>% unique), "Set3")
  df_plot <- geno_df 
  
  p1 <- ggplot(df_plot %>% group_by(!!as.name(col_to_group),stage , drug) %>% mutate(N = length(unique(key)), group = paste0(!!as.name(col_to_group), " (", N ,")") ) %>% ungroup() %>% 
                 select(-key, -weights, -!!as.name(col_to_group)) %>% unique() %>% filter(N > !!min_n) , aes(x =group , y = avg_prob, fill = archs )) + 
    geom_col() + theme_bw() + labs(xlab = "%") + facet_wrap(stage ~ drug, scales = "free") + RotatedAxis()   + scale_fill_manual("Archs", values  = palette_archs)
  
  p2 <- ggplot(df_plot %>% group_by(!!as.name(col_to_group), stage) %>% mutate(N = length(unique(key)), group = paste0(!!as.name(col_to_group), " (", N ,")") ) %>% ungroup() %>% 
                 select(-key, -weights, -!!as.name(col_to_group)) %>% unique()  %>% filter(N > !!min_n), aes(x = group, y = avg_prob, fill = archs )) + 
    geom_col() + theme_bw() + labs(xlab = "%") + facet_wrap(stage ~ ., scales = "free") + RotatedAxis()   + scale_fill_manual("Archs", values  = palette_archs)
  p3 <- ggplot(df_plot %>% group_by(!!as.name(col_to_group)) %>% mutate(N = length(unique(key)), group = paste0(!!as.name(col_to_group), " (", N ,")") ) %>% ungroup() %>% 
                 select(-key, -weights, -!!as.name(col_to_group), drug) %>% unique()  %>% filter(N > !!min_n) , aes(x =group, y = avg_prob, fill = archs )) + 
    geom_col() + theme_bw() + labs(xlab = "%") + facet_wrap(drug ~ ., scales = "free") + RotatedAxis()   + scale_fill_manual("Archs", values  = palette_archs)
  p4 <- ggplot(df_plot %>% group_by(!!as.name(col_to_group)) %>% mutate(N = length(unique(key)), group = paste0(!!as.name(col_to_group), " (", N ,")") ) %>% ungroup() %>% 
                 select(-key, -weights, -!!as.name(col_to_group)) %>% unique()  %>% filter(N > !!min_n), aes(x =group, y = avg_prob, fill = archs )) + 
    geom_col() + theme_bw() + labs(xlab = "%")  + RotatedAxis()   + scale_fill_manual("Archs", values  = palette_archs)
  # p2 <-  ggplot(geno_df %>% select(-key, -weights) %>% unique(), aes(x = archs, y = avg_prob, fill = !!as.name(col_to_group) )) + 
  #   geom_col() + theme_bw() + labs(xlab = "N") + facet_wrap(stage ~ drug, scales = "free") + RotatedAxis()   + scale_fill_manual("Clusters", values  = palette)
  
  return(list("by_all" = p1,"by_stage" = p2, "by_drug" = p3, "by_none" = p4))
}




from_triple_to_ridgeplot <- function(df_barplot,group_by_col,  bc44_inp, drug_inp,title_add,
                                     title_add_col ,title_colour = "black", label = c("Parental", "Parental", "Under drug", "Re-growth")) {
  
  library(ggridges)
  library(ggalluvial)
  library(tidyverse)
  library(ggtext)
  library(RColorBrewer)
  
  # Auto-generate archetype colors
  generate_colors <- function(n) {
    brewer_palette <- brewer.pal(n, "Dark2") # Get the max number of colors (11) from a diverging palette
    #color_func <- colorRampPalette(brewer_palette) # Interpolate these colors
    return(brewer_palette) # Generate n colors
  }
  archetype_colors <- generate_colors(length( unique(df_barplot$archs))) 
  
  names(archetype_colors) <- unique(df_barplot$archs)
  
  df_barplot$stage <- factor(df_barplot$stage, levels = c("parental", "under_drug", "expansion"))
  
  
  title_text <- glue::glue("<span style='color:{title_colour};'>{title_add_col} </span> ")
  
  title_add_col_axis <- gsub(" barcode", "", title_add_col)
  axis_to_add <- glue::glue("<span style='color:{title_colour};'>({title_add_col_axis}) </span> ")
  
  ggplot(df_barplot %>% filter((!!as.name(group_by_col) == !!bc44_inp) |
                                 ((!!as.name(group_by_col) == "other") & (stage == "parental")), 
                               drug == !!drug_inp  | stage == "parental") %>% 
           ungroup() %>% 
           dplyr::mutate(stage = if_else((stage == "parental") & (!!as.name(group_by_col) == "other"),
                                         "parental_other", stage)) %>% 
           select(avg_prob, archs, stage) %>% unique() %>% 
           mutate(stage = factor(stage, levels = c("parental_other","parental", "under_drug", "expansion"))),
         aes(y = avg_prob, x = stage)) +
    geom_flow(aes(alluvium = archs), alpha= .9, 
              lty = 2, fill = "white", color = "black",
              curve_type = "linear", 
              width = .5) +
    geom_col(aes(fill = archs), width = .5, color = "black") +
    scale_fill_manual("", 
                      values= archetype_colors)+
    scale_y_continuous(expand = c(0,0)) +
    cowplot::theme_minimal_hgrid(font_size = 15) + xlab("") + Seurat::RotatedAxis() + 
    ylab("AA weight") +
    labs(title = paste0('Avg. composition ', title_text, title_add)) + 
    scale_x_discrete(labels = c(glue::glue("{label[1]} (others)"),
                                glue::glue("{label[2]} {axis_to_add}") ,
                                glue::glue("{label[3]} {axis_to_add}") , 
                                glue::glue("{label[4]} {axis_to_add}"))
    ) +
    theme( plot.title = element_markdown(), text = element_text(family = "Arial"), 
           axis.text.x = element_markdown() )
  
  
}



from_double_to_ridgeplot <- function(df_barplot,group_by_col,  bc44_inp, drug_inp,title_add,
                                     title_add_col ,title_colour = "black", label = c("Parental", "Parental", "Re-growth")) {
  
  library(ggridges)
  library(ggalluvial)
  library(tidyverse)
  library(ggtext)
  library(RColorBrewer)
  
  # Auto-generate archetype colors
  generate_colors <- function(n) {
    brewer_palette <- brewer.pal(n, "Set2") # Get the max number of colors (11) from a diverging palette
    #color_func <- colorRampPalette(brewer_palette) # Interpolate these colors
    return(brewer_palette) # Generate n colors
  }
  archetype_colors <- generate_colors(length( unique(df_barplot$archs))) 
  
  names(archetype_colors) <- unique(df_barplot$archs)
  
  df_barplot$stage <- factor(df_barplot$stage, levels = c("parental",  "expansion"))
  
  
  title_text <- glue::glue("<span style='color:{title_colour};'>{title_add_col} </span> ")
  
  title_add_col_axis <- gsub(" barcode", "", title_add_col)
  axis_to_add <- glue::glue("<span style='color:{title_colour};'>({title_add_col_axis}) </span> ")
  
  ggplot(df_barplot %>% filter((!!as.name(group_by_col) == !!bc44_inp) |
                                 ((!!as.name(group_by_col) == "other") & (stage == "parental")), 
                               drug == !!drug_inp  | stage == "parental") %>% 
           ungroup() %>% 
           dplyr::mutate(stage = if_else((stage == "parental") & (!!as.name(group_by_col) == "other"),
                                         "parental_other", stage)) %>% 
           select(avg_prob, archs, stage) %>% unique() %>% 
           mutate(stage = factor(stage, levels = c("parental_other","parental", "expansion"))),
         aes(y = avg_prob, x = stage)) +
    geom_flow(aes(alluvium = archs), alpha= .9, 
              lty = 2, fill = "white", color = "black",
              curve_type = "linear", 
              width = .5) +
    geom_col(aes(fill = archs), width = .5, color = "black") +
    scale_fill_manual("", 
                      values= archetype_colors)+
    scale_y_continuous(expand = c(0,0)) +
    cowplot::theme_minimal_hgrid(font_size = 15) + xlab("") + Seurat::RotatedAxis() + 
    ylab("AA weight") +
    labs(title = paste0('Avg. composition ', title_text, title_add)) + 
    scale_x_discrete(labels = c(glue::glue("{label[1]} (others)"),
                                glue::glue("{label[2]} {axis_to_add}"), 
                                glue::glue("{label[3]} {axis_to_add}"))
    ) +
    theme( plot.title = element_markdown(), text = element_text(family = "Arial"), 
           axis.text.x = element_markdown() )
  
  
}


add_rescaled_other_bc <- function(df, bc, new_name = "other") {
  library(dplyr)
  df %>% dplyr::mutate(bc_new = if_else(bc44 == !!bc,!!bc, !!new_name)) %>% 
    dplyr::group_by(across(c(-avg_prob, -bc44))) %>% 
    dplyr::summarize(avg_prob = mean(avg_prob))
}

add_scPubr_theme <- function(p){
  library(ggplot2)
  p <- p &
    ggplot2::theme_minimal(base_size = 12) &
    ggplot2::theme(
      plot.title.position = "plot",
      plot.caption.position = "plot",
      text = ggplot2::element_text(family = "Arial"),
      legend.justification = "center",
      legend.position = "bottom",
      panel.grid = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0),
      plot.background = ggplot2::element_rect(fill = "white", color = "white"),
      panel.background = ggplot2::element_rect(fill = "white", color = "white"),
      legend.background = ggplot2::element_rect(fill = "white", color = "white"),
      axis.text.x = element_blank(),
      axis.text.y = element_blank()) & ylab("") & xlab("")
  p
}

