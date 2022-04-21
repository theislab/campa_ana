# clear the workspace
rm(list=ls())

# load required packages
library(tidyverse)

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

# directory to hold results of mixed models
model_dir <- file.path(campa_ana$constants$SOURCE_DIR,"mixed_model_results")

# create directory to hold plots
plot_dir <- file.path(campa_ana$constants$SOURCE_DIR,"figures","dot_plots")
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
}

# read fold changes calculated previously
intensity_fold_changes <- read_csv(file = file.path(model_dir,"184A1_intensity_fold_changes.csv"))
size_fold_changes <- read_csv(file = file.path(model_dir,"184A1_size_fold_changes.csv"))

# make bubble plot
make_bubble_plot <- function(fold_changes,
                             plot_var,
                             row_var,
                             col_var,
                             row_clustering=NULL,
                             col_clustering=NULL,
                             color_limits=c(-2,2)) {
  
  plot_var = enquo(plot_var)
  row_var = enquo(row_var)
  col_var = enquo(col_var)
  
  if (!is.null(row_clustering)) {
    fold_changes <- mutate(fold_changes,!!row_var:=factor(!!row_var,levels=row_clustering))
  } 
  if (!is.null(col_clustering)) {
    fold_changes <- mutate(fold_changes,!!col_var:=factor(col_var,levels=col_clustering))
  }
  
  bubble_plot <- fold_changes %>%
    mutate(
      p.value.adj.BY = p.adjust(p.value, "BY"),
      signif_code = case_when(
        p.value.adj.BY > 0.05 ~ 0,
        p.value.adj.BY > 0.01 ~ 1,
        TRUE ~ 2)) %>%
    
    ggplot(aes(y=!!col_var,
               x=!!row_var,
               col=!!plot_var,
               size=signif_code)) +
    geom_point() +
    scale_radius(name = "",
                 breaks = c(0,1,2),
                 labels = c("p > 0.05","p < 0.05","p < 0.01"),
                 range = c(0.5,2.5),
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

intensity_fold_changes %>%
  #filter(cluster=="Nuclear speckles") %>%
  filter(treatment=="Meayamycin (12.5h) - Control") %>%
  filter(comparison=="unnormalised") %>%
  mutate(
    p.value.adj.BY = p.adjust(p.value, "BY")) %>% 
  #filter(channel=="01_PABPC1") %>%
  #filter(cluster=="Nuclear speckles - all")
  make_bubble_plot(plot_var = log2_fold_change,
                   row_var = channel,
                   col_var = cluster)
