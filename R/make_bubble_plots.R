# clear the workspace
rm(list=ls())

# load required packages
library(tidyverse)
library(Cairo)
library(patchwork)

# setup script-specific parameters
cell_type <- "184A1"
experiment_name <- "VAE_all/CondVAE_pert-CC"

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

# read fold changes calculated previously 
intensity_fold_changes_direct_norm <- map_dfr(
  list.files(model_dir) %>% 
    grepl(pattern = "intensity_fold_changes_direct_norm") %>%
    list.files(model_dir)[.],
  ~read_csv(file = file.path(model_dir,.))
)

intensity_fold_changes_full <- map_dfr(
  list.files(model_dir) %>% 
    grepl(pattern = "intensity_fold_changes_(?!(direct_norm|DMSO))",x = .,perl=T) %>%
    list.files(model_dir)[.],
  ~read_csv(file = file.path(model_dir,.))
)

size_fold_changes_direct_norm <- map_dfr(
  list.files(model_dir) %>% 
    grepl(pattern = "size_fold_changes_direct_norm") %>%
    list.files(model_dir)[.],
  ~read_csv(file = file.path(model_dir,.))
)

size_fold_changes_full <- map_dfr(
  list.files(model_dir) %>% 
    grepl(pattern = "size_fold_changes_(?!(direct_norm|DMSO))",x = .,perl=T) %>%
    list.files(model_dir)[.],
  ~read_csv(file = file.path(model_dir,.))
)

# Compare the two model fitting approaches ----

# compare unnormalised results
unnormalised_comparison <- left_join(
  intensity_fold_changes_full %>%
  filter(comparison=="unnormalised") %>% 
  mutate(treatment = str_remove(treatment," - Unperturbed")),
  intensity_fold_changes_direct_norm %>%
    filter(comparison=="unnormalised"),
  by=c("treatment","comparison","channel","cluster"),
  suffix=c("_full","_direct_norm")) %>%
  filter(cluster!="all")

compare_fold_changes <- unnormalised_comparison %>%
  ggplot(aes(x=fold_change_direct_norm,y=fold_change_full,col=treatment)) + 
  geom_abline(size=0.2) + 
  geom_point(size=0.2) + 
  scale_x_log10() + 
  scale_y_log10() +
  geom_vline(xintercept = 1,size=0.2) +
  geom_hline(yintercept = 1,size=0.2) +
  theme_bw(base_size = 8) +
  theme(panel.grid = element_blank(),
        legend.title=element_blank())

compare_p_values <- unnormalised_comparison %>%
  ggplot(aes(x=p.value_direct_norm,y=p.value_full,col=treatment)) + 
  geom_abline(size=0.2) + 
  geom_point(size=0.2) + 
  scale_x_log10() + 
  scale_y_log10() +
  geom_vline(xintercept = 0.05,size=0.2) +
  geom_hline(yintercept = 0.05,size=0.2) +
  theme_bw(base_size = 8) +
  theme(panel.grid = element_blank(),
        legend.title=element_blank())


CairoPNG(filename=file.path(plot_dir,"direct_full_comparison_unnormalised_fold_changes.png"),height=7,width=10,units="cm",res=600)
print(ggExtra::ggMarginal(compare_fold_changes,type = "histogram",margins = "x",bins = 100))
dev.off()

CairoPNG(filename=file.path(plot_dir,"direct_full_comparison_unnormalised_p_values.png"),height=7,width=10,units="cm",res=600)
print(ggExtra::ggMarginal(compare_p_values,type = "histogram",margins = "x",bins = 100))
dev.off()

# compare normalised results
normalised_comparison <- left_join(
  intensity_fold_changes_full %>%
    filter(comparison=="relative_to_all") %>%
    mutate(treatment = str_remove(treatment," - Unperturbed"),
           cluster = str_remove(cluster," - all")),
  intensity_fold_changes_direct_norm %>%
    filter(comparison=="relative_to_all"),
  by=c("treatment","comparison","channel","cluster"),
  suffix=c("_full","_direct_norm"))

compare_fold_changes <- normalised_comparison %>%
  ggplot(aes(x=fold_change_direct_norm,y=fold_change_full,col=treatment)) + 
  geom_abline(size=0.2) + 
  geom_point(size=0.2) + 
  scale_x_log10() + 
  scale_y_log10() +
  geom_vline(xintercept = 1,size=0.2) +
  geom_hline(yintercept = 1,size=0.2) +
  theme_bw(base_size = 8) +
  theme(panel.grid = element_blank(),
        legend.title=element_blank())

compare_p_values <- normalised_comparison %>%
  ggplot(aes(x=p.value_direct_norm,y=p.value_full,col=treatment)) + 
  geom_abline(size=0.2) + 
  geom_point(size=0.2) + 
  scale_x_log10() + 
  scale_y_log10() +
  geom_vline(xintercept = 0.05,size=0.2) +
  geom_hline(yintercept = 0.05,size=0.2) +
  theme_bw(base_size = 8) +
  theme(panel.grid = element_blank(),
        legend.title=element_blank())


CairoPNG(filename=file.path(plot_dir,"direct_full_comparison_normalised_fold_changes.png"),height=7,width=10,units="cm",res=600)
print(ggExtra::ggMarginal(compare_fold_changes,type = "histogram",margins = "x",bins = 100))
dev.off()

CairoPNG(filename=file.path(plot_dir,"direct_full_comparison_normalised_p_values.png"),height=7,width=10,units="cm",res=600)
print(ggExtra::ggMarginal(compare_p_values,type = "histogram",margins = "x",bins = 100))
dev.off()


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


# note that "all" is computed every time there is an unnormalised pairwise comparison, just keep the maximum p.value (conservative)

fold_changes_to_plot <- intensity_fold_changes_full %>%
  filter(cluster=="all") %>%
  group_by(cluster,channel,treatment) %>%
  filter(p.value==max(p.value)) %>%
  ungroup() %>%
  bind_rows(filter(intensity_fold_changes_full,cluster!="all")) %>%
  mutate(treatment = str_remove(treatment," - Unperturbed"),
         cluster = str_remove(cluster," - all")) %>%
  left_join(rename(rename_channels_for_plotting_dict,channel_name=new),
            by=c("channel"="old")) %>%
  mutate(cluster=if_else(cluster=="all","Whole nucleus",cluster))

size_fold_changes_to_plot <- size_fold_changes_full %>%
  filter(cluster=="all") %>%
  group_by(cluster,channel,treatment) %>%
  filter(p.value==max(p.value)) %>%
  ungroup() %>%
  bind_rows(filter(size_fold_changes_full,cluster!="all")) %>%
  mutate(treatment = str_remove(treatment," - Unperturbed"),
         cluster = str_remove(cluster," - all")) %>%
  left_join(rename(rename_channels_for_plotting_dict,channel_name=new),
            by=c("channel"="old")) %>%
  mutate(cluster=if_else(cluster=="all","Whole nucleus",cluster))

# 34 channels x (6 clusters + overall) x 6 treatments = 1428
fold_changes_to_plot %>%
  filter(comparison=="unnormalised") %>%
  nrow()

# (6 clusters + overall) x 6 treatments = 42
size_fold_changes_to_plot %>%
  filter(comparison=="unnormalised") %>%
  nrow()

# 34 channels x 6 clusters x 6 treatments = 1224
fold_changes_to_plot %>%
  filter(comparison=="relative_to_all") %>%
  nrow()

# 6 clusters x 6 treatments = 1224
size_fold_changes_to_plot %>%
  filter(comparison=="relative_to_all") %>%
  nrow()

# Loop over all treatments
all_treatments <- unique(fold_changes_to_plot$treatment)
for (current_treatment in all_treatments) {
  
  # Plot unnormalised
  to_plot <- fold_changes_to_plot %>%
    filter(treatment==current_treatment) %>%
    filter(comparison=="unnormalised")
  cols <- hclust_cols(to_plot)
  rows <- hclust_rows(to_plot)
  sizes_to_plot <- size_fold_changes_to_plot %>%
    filter(treatment==current_treatment) %>%
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
  
  combined_plot <- (size_bubble_plot  +
                      theme(axis.title.x = element_blank(),
                            axis.title.y = element_blank(),
                            legend.position = "none") ) + 
    ( bubble_plot +
        xlab("Channel") + 
        theme(axis.text.y=element_blank(),
              axis.ticks.y=element_blank()) ) + 
    plot_layout(widths=c(1,22))
  combined_plot
  
  ggsave(plot = combined_plot,filename = file.path(plot_dir,paste0("unnormalised_dot_plot_",current_treatment,".pdf")),width=10,height=5,units="cm")
  ggsave_cairo(plot = combined_plot,filename = file.path(plot_dir,paste0("unnormalised_dot_plot_",current_treatment,".png")),width=10,height=5,units="cm",dpi=600)
  
  
  # Plot normalised
  to_plot <- fold_changes_to_plot %>%
    filter(treatment==current_treatment) %>%
    filter(comparison=="relative_to_all")
  sizes_to_plot <- size_fold_changes_to_plot %>%
    filter(treatment==current_treatment) %>%
    filter(comparison=="relative_to_all") %>%
    mutate(axis_label="Fraction of\npixels")
  
  bubble_plot <- make_bubble_plot(fold_changes = to_plot,
                                  plot_var = log2_fold_change,
                                  col_var = channel_name,
                                  row_var = cluster,
                                  col_order=cols,
                                  row_order=rows,
                                  color_limits=c(-1,1))
  
  size_bubble_plot <- make_bubble_plot(fold_changes = sizes_to_plot,
                                       plot_var = log2_fold_change,
                                       col_var = axis_label,
                                       row_var = cluster,
                                       row_order=rows,
                                       color_limits=c(-1,1))
  
  combined_plot <- (size_bubble_plot  +
                      theme(axis.title.x = element_blank(),
                            axis.title.y = element_blank(),
                            legend.position = "none") ) + 
    ( bubble_plot +
        xlab("Channel") + 
        theme(axis.text.y=element_blank(),
              axis.ticks.y=element_blank()) ) + 
    plot_layout(widths=c(1,22))
  combined_plot
  
  ggsave(plot = combined_plot,filename = file.path(plot_dir,paste0("normalised_dot_plot_",current_treatment,".pdf")),width=10,height=5,units="cm")
  ggsave_cairo(plot = combined_plot,filename = file.path(plot_dir,paste0("normalised_dot_plot_",current_treatment,".png")),width=10,height=5,units="cm",dpi=600)
  
}
