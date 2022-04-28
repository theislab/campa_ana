# clear the workspace
rm(list=ls())

# load required packages
library(tidyverse)
library(Cairo)
library(patchwork)

# setup script-specific parameters
cell_type <- "HeLa"
experiment_name <- "VAE_SBF2/CondVAE_siRNA-CC"

# setup python environment and experiment-specific parameters
library(reticulate)
reticulate::use_condaenv("pelkmans-3.9")
campa <- import("campa")
campa_ana <- import("campa_ana")
np <- import("numpy") 

# load required R functions from the campa_ana package
source(file.path(campa_ana$constants$SOURCE_DIR,"R","setup_paths.R"))
source(file.path(campa_ana$constants$SOURCE_DIR,"R","io.R"))
source(file.path(campa_ana$constants$SOURCE_DIR,"R","mixed_models.R"))
source(file.path(campa_ana$constants$SOURCE_DIR,"R","channels_rename_info.R"))

# directory to hold results of mixed models
model_dir <- file.path(campa_ana$constants$SOURCE_DIR,"mixed_model_results")

# create directory to hold plots
plot_dir <- file.path(campa_ana$constants$SOURCE_DIR,"figures","dot_plots")
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
}

intensity_fold_changes_full <- read_csv(file = file.path(model_dir,"HeLa_intensity_fold_changes_SBF2.csv"))
size_fold_changes_full <- read_csv(file = file.path(model_dir,"HeLa_size_fold_changes_SBF2.csv"))

nuclear_CSLs <- c(
  "Cajal bodies",             
  "Nuclear periphery",        
  "Nuclear speckles",         
  "Nucleolus",                
  "Nucleoplasm",
  "PML bodies")

# Make bubble plots ----

make_bubble_plot <- function(fold_changes,
                             plot_var,
                             row_var,
                             col_var,
                             row_order=NULL,
                             col_order=NULL,
                             color_limits=c(-2,2)) {
  
  plot_var = enquo(plot_var)
  row_var = enquo(row_var)
  col_var = enquo(col_var)
  
  if (!is.null(col_order)) {
    fold_changes <- mutate(fold_changes,!!col_var:=factor(!!col_var,levels=col_order))
  } 
  if (!is.null(row_order)) {
    fold_changes <- mutate(fold_changes,!!row_var:=factor(!!row_var,levels=row_order))
  }
  
  bubble_plot <- fold_changes %>%
    mutate(
      p.value.adj.BY = p.adjust(p.value, "BY"),
      signif_code = case_when(
        p.value.adj.BY > 0.05 ~ 0,
        p.value.adj.BY > 0.01 ~ 1,
        TRUE ~ 2)) %>%
    
    ggplot(aes(x=!!col_var,
               y=!!row_var,
               col=!!plot_var,
               size=signif_code)) +
    geom_point() +
    scale_radius(name = "",
                 breaks = c(0,1,2),
                 labels = c("p > 0.05","p < 0.05","p < 0.01"),
                 range = c(0.75,1.75),
                 limits = c(0,2),
                 guide = guide_legend(order=2)) +
    scale_x_discrete(position = "bottom") +
    scale_y_discrete(position = "left") +
    scale_colour_gradient2(name= expression(log[2](`fold-change`)),
                           breaks=c(color_limits[1],0,color_limits[2]),
                           low = "blue",mid="white",high = "red",midpoint = 0,
                           limits=color_limits,oob=scales::squish,
                           guide = guide_colourbar(order=1)) +
    theme_bw(base_size = 7) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle=90,hjust=1,vjust=0.5),
          axis.text.y = element_text(),
          axis.title.y = element_blank(),
          legend.key.size = unit(0.15,"cm"),
          legend.title = element_text(size=6),
          legend.box = "horizontal",
          legend.position="bottom",
          legend.box.margin = margin(0,0,0,0,unit="cm"))
  bubble_plot
}

hclust_rows <- function(df) {
  row_clustering <- df %>%
    select(cluster,channel_name,fold_change) %>%
    pivot_wider(names_from = channel_name,values_from=fold_change) %>%
    column_to_rownames("cluster") %>%
    dist() %>%
    hclust(method = "complete")
  row_order <- order.dendrogram(as.dendrogram(row_clustering))
  row_order_names <- c(setdiff(row_clustering$labels[row_order],"Whole nucleus"),"Whole nucleus")
  return(row_order_names)
}

hclust_cols <- function(df) {
  col_clustering <- df %>%
    select(cluster,channel_name,fold_change) %>%
    pivot_wider(names_from = channel_name,values_from=fold_change) %>%
    column_to_rownames("cluster") %>%
    t() %>%
    dist() %>%
    hclust(method = "complete")
  col_order <- order.dendrogram(as.dendrogram(col_clustering))
  col_order_names <- col_clustering$labels[col_order]
  return(col_order_names)
}

# note that "Nucleus (combined)","Cytoplasm (combined)" are computed every time
# there is an unnormalised pairwise comparison, just keep the maximum p.value
# (conservative)

fold_changes_to_plot <- intensity_fold_changes_full %>%
  filter(cluster %in% c("Nucleus (combined)","Cytoplasm (combined)")) %>%
  group_by(cluster,channel,siRNA) %>%
  filter(p.value==max(p.value)) %>%
  ungroup() %>%
  bind_rows(filter(intensity_fold_changes_full,!(cluster %in% c("Nucleus (combined)","Cytoplasm (combined)")))) %>%
  mutate(siRNA = str_remove(siRNA," - scrambled"),
         cluster = str_remove(cluster," - Nucleus \\(combined\\)"),
         cluster = str_remove(cluster," - Cytoplasm \\(combined\\)")) %>%
  left_join(rename(rename_channels_for_plotting_dict,channel_name=new),
            by=c("channel"="old")) %>%
  mutate(cluster = if_else(substr(cluster,1,1)=="(", 
                           str_remove_all(cluster,"^\\(|\\)$"),
                           cluster))

size_fold_changes_to_plot <- size_fold_changes_full %>%
  filter(cluster %in% c("Nucleus (combined)","Cytoplasm (combined)")) %>%
  group_by(cluster,channel,siRNA) %>%
  filter(p.value==max(p.value)) %>%
  ungroup() %>%
  bind_rows(filter(size_fold_changes_full,!(cluster %in% c("Nucleus (combined)","Cytoplasm (combined)")))) %>%
  mutate(siRNA = str_remove(siRNA," - scrambled"),
         cluster = str_remove(cluster," - Nucleus \\(combined\\)"),
         cluster = str_remove(cluster," - Cytoplasm \\(combined\\)")) %>%
  left_join(rename(rename_channels_for_plotting_dict,channel_name=new),
            by=c("channel"="old")) %>%
  mutate(cluster = if_else(substr(cluster,1,1)=="(", 
                           str_remove_all(cluster,"^\\(|\\)$"),
                           cluster))

fold_changes_to_plot %>%
  group_by(siRNA,cluster,comparison) %>% count()

# 43 channels x (15 clusters + nucleus + cytoplasm) = 731
fold_changes_to_plot %>%
  filter(comparison=="unnormalised") %>%
  nrow()

# (15 clusters + nucleus + cytoplasm) = 17
size_fold_changes_to_plot %>%
  filter(comparison=="unnormalised") %>%
  nrow()

# 43 channels x 6 clusters = 258
fold_changes_to_plot %>%
  filter(comparison=="relative_to_Nucleus (combined)") %>%
  nrow()

# 6 clusters = 6
size_fold_changes_to_plot %>%
  filter(comparison=="relative_to_Nucleus (combined)") %>%
  nrow()

# 43 channels x 9 clusters = 387
fold_changes_to_plot %>%
  filter(comparison=="relative_to_Cytoplasm (combined)") %>%
  nrow()

# 9 clusters = 9
size_fold_changes_to_plot %>%
  filter(comparison=="relative_to_Cytoplasm (combined)") %>%
  nrow()

current_siRNA <- "SBF2"

# Plot unnormalised
to_plot <- fold_changes_to_plot %>%
  filter(siRNA==current_siRNA) %>%
  filter(comparison=="unnormalised")
cols <- hclust_cols(to_plot)
rows <- hclust_rows(to_plot)
sizes_to_plot <- size_fold_changes_to_plot %>%
  filter(siRNA==current_siRNA) %>%
  filter(comparison=="unnormalised") %>%
  mutate(axis_label="Number of\npixels")

bubble_plot <- make_bubble_plot(fold_changes = to_plot,
                                plot_var = log2_fold_change,
                                col_var = channel_name,
                                row_var = cluster,
                                col_order=cols,
                                row_order=rows,
                                color_limits=c(-2,2))

size_bubble_plot <- make_bubble_plot(fold_changes = sizes_to_plot,
                                     plot_var = log2_fold_change,
                                     col_var = axis_label,
                                     row_var = cluster,
                                     row_order=rows,
                                     color_limits=c(-2,2))

combined_plot <- (
  size_bubble_plot +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none")
) + (
  bubble_plot +
    xlab("Channel") + 
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) 
) + 
  plot_layout(widths=c(1,30))
combined_plot

ggsave(plot = combined_plot,filename = file.path(plot_dir,paste0("unnormalised_dot_plot_",current_siRNA,".pdf")),width=16,height=8,units="cm")
ggsave_cairo(plot = combined_plot,filename = file.path(plot_dir,paste0("unnormalised_dot_plot_",current_siRNA,".png")),width=16,height=8,units="cm",dpi=600)

# Plot normalised
to_plot_nucleus <- fold_changes_to_plot %>%
  filter(siRNA==current_siRNA) %>%
  filter(comparison=="relative_to_Nucleus (combined)")
to_plot_cytoplasm <- fold_changes_to_plot %>%
  filter(siRNA==current_siRNA) %>%
  filter(comparison=="relative_to_Cytoplasm (combined)")
sizes_to_plot_nucleus <- size_fold_changes_to_plot %>%
  filter(siRNA==current_siRNA) %>%
  filter(comparison=="relative_to_Nucleus (combined)") %>%
  mutate(axis_label="Fraction of\npixels")
sizes_to_plot_cytoplasm <- size_fold_changes_to_plot %>%
  filter(siRNA==current_siRNA) %>%
  filter(comparison=="relative_to_Cytoplasm (combined)") %>%
  mutate(axis_label="Fraction of\npixels")

nucleus_bubble_plot <- make_bubble_plot(fold_changes = to_plot_nucleus,
                                plot_var = log2_fold_change,
                                col_var = channel_name,
                                row_var = cluster,
                                col_order=cols,
                                row_order=rows,
                                color_limits=c(-1,1))

nucleus_size_bubble_plot <- make_bubble_plot(fold_changes = sizes_to_plot_nucleus,
                                     plot_var = log2_fold_change,
                                     col_var = axis_label,
                                     row_var = cluster,
                                     row_order=rows,
                                     color_limits=c(-1,1))

cytoplasm_bubble_plot <- make_bubble_plot(fold_changes = to_plot_cytoplasm,
                                        plot_var = log2_fold_change,
                                        col_var = channel_name,
                                        row_var = cluster,
                                        col_order=cols,
                                        row_order=rows,
                                        color_limits=c(-1,1))

cytoplasm_size_bubble_plot <- make_bubble_plot(fold_changes = sizes_to_plot_cytoplasm,
                                             plot_var = log2_fold_change,
                                             col_var = axis_label,
                                             row_var = cluster,
                                             row_order=rows,
                                             color_limits=c(-1,1))

combined_plot <- (nucleus_size_bubble_plot  +
                    theme(axis.title.x = element_blank(),
                          axis.title.y = element_blank(),
                          axis.text.x=element_blank(),
                          axis.ticks.x=element_blank(),
                          legend.position = "none") ) + 
  ( nucleus_bubble_plot  +
      theme(axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.x=element_blank()) ) +
  ( cytoplasm_size_bubble_plot  +
     theme(axis.title.x = element_blank(),
           axis.title.y = element_blank(),
           legend.position = "none") ) + 
  ( cytoplasm_bubble_plot +
      xlab("Channel") + 
      theme(axis.text.y=element_blank(),
            axis.ticks.y=element_blank()) ) +
  plot_layout(widths=c(1,30),heights=c(6,9),ncol = 2,nrow=2,guides="collect") & theme(legend.position = "bottom")
combined_plot

ggsave(plot = combined_plot,filename = file.path(plot_dir,paste0("normalised_dot_plot_",current_siRNA,".pdf")),width=16,height=9,units="cm")
ggsave_cairo(plot = combined_plot,filename = file.path(plot_dir,paste0("normalised_dot_plot_",current_siRNA,".png")),width=16,height=9,units="cm",dpi=600)

# make just the size bubble plots for the main figure

row_order_nucleus <- c("Nucleus (combined)",rows[rows %in% nuclear_CSLs])
row_order_cytoplasm <- c(rows[!rows %in% row_order_nucleus])

sizes_to_plot_nucleus_un_norm <- size_fold_changes_to_plot %>%
  filter(cluster %in% row_order_nucleus) %>%
  filter(comparison=="unnormalised") %>%
  select(cluster,p.value,log2_fold_change) %>%
  mutate(axis_label="Number of\npixels")

sizes_to_plot_cytoplasm_un_norm <- size_fold_changes_to_plot %>%
  filter(cluster %in% row_order_cytoplasm) %>%
  filter(comparison=="unnormalised") %>%
  select(cluster,p.value,log2_fold_change) %>%
  mutate(axis_label="Number of\npixels") 

sizes_to_plot_nucleus_norm <- size_fold_changes_to_plot %>%
  filter(cluster %in% row_order_nucleus) %>%
  filter(comparison=="relative_to_Nucleus (combined)") %>%
  select(cluster,p.value,log2_fold_change) %>%
  bind_rows(tibble(cluster="Nucleus (combined)",log2_fold_change=0)) %>%
  mutate(axis_label="Fraction of\npixels")

sizes_to_plot_cytoplasm_norm <- size_fold_changes_to_plot %>%
  filter(cluster %in% row_order_cytoplasm) %>%
  filter(comparison=="relative_to_Cytoplasm (combined)") %>%
  select(cluster,p.value,log2_fold_change) %>%
  bind_rows(tibble(cluster="Cytoplasm (combined)",log2_fold_change=0)) %>%
  mutate(axis_label="Fraction of\npixels")


nucleus_sizes_bubble_plot <- make_bubble_plot(
  fold_changes = sizes_to_plot_nucleus_un_norm,
  plot_var = log2_fold_change,
  col_var = axis_label,
  row_var = cluster,
  row_order=row_order_nucleus,
  color_limits=c(-2,2))

nucleus_sizes_bubble_plot_norm <- make_bubble_plot(
  fold_changes = sizes_to_plot_nucleus_norm,
  plot_var = log2_fold_change,
  col_var = axis_label,
  row_var = cluster,
  row_order=row_order_nucleus,
  color_limits=c(-2,2))

cytoplasm_sizes_bubble_plot <- make_bubble_plot(
  fold_changes = sizes_to_plot_cytoplasm_un_norm,
  plot_var = log2_fold_change,
  col_var = axis_label,
  row_var = cluster,
  row_order=row_order_cytoplasm,
  color_limits=c(-2,2))

cytoplasm_sizes_bubble_plot_norm <- make_bubble_plot(
  fold_changes = sizes_to_plot_cytoplasm_norm,
  plot_var = log2_fold_change,
  col_var = axis_label,
  row_var = cluster,
  row_order=row_order_cytoplasm,
  color_limits=c(-2,2))

combined_plot <- nucleus_sizes_bubble_plot + theme(axis.text.x=element_blank(),
                                                   axis.title.x=element_blank(),
                                                   axis.ticks.x=element_blank()) +
  nucleus_sizes_bubble_plot_norm + theme(axis.text.x=element_blank(),
                                         axis.title.x=element_blank(),
                                         axis.ticks.x=element_blank(),
                                         axis.text.y=element_blank(),
                                         axis.title.y=element_blank(),
                                         axis.ticks.y=element_blank()) +
  cytoplasm_sizes_bubble_plot + theme(axis.title.x=element_blank()) +
  cytoplasm_sizes_bubble_plot_norm + theme(axis.title.x=element_blank(),
                                           axis.text.y=element_blank(),
                                           axis.title.y=element_blank(),
                                           axis.ticks.y=element_blank()) +
  plot_layout(widths=c(1,1),heights=c(6,9),ncol = 2,nrow=2,guides="collect") & theme(legend.position = "none")

combined_plot

ggsave(plot = combined_plot,filename = file.path(plot_dir,paste0("normalised_dot_plot_sizes_only_",current_siRNA,".pdf")),width=4.2,height=6.5,units="cm")
ggsave_cairo(plot = combined_plot,filename = file.path(plot_dir,paste0("normalised_dot_plot_sizes_only_",current_siRNA,".png")),width=4.2,height=6.5,units="cm",dpi=600)

ggpubr::as_ggplot(
  ggpubr::get_legend(combined_plot & theme(legend.position = "right"))
)
ggsave(filename = file.path(plot_dir,paste0("normalised_dot_plot_sizes_only_",current_siRNA,"_legend.pdf")),width=4.2,height=3,units="cm")
ggsave_cairo(filename = file.path(plot_dir,paste0("normalised_dot_plot_sizes_only_",current_siRNA,"_legend.png")),width=4.2,height=3,units="cm",dpi=600)

