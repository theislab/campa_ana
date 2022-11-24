# clear the workspace
rm(list=ls())

# load required packages
library(tidyverse)
library(ggtree)

# setup script-specific parameters
cell_type <- "184A1"
experiment_name <- "VAE_all/CondVAE_pert-CC"

# setup python environment and experiment-specific parameters
library(reticulate)
reticulate::use_condaenv("pelkmans-3.9")
#campa <- import("campa")
campa_ana <- import("campa_ana")

# load required R functions from the campa_ana package
source(file.path(campa_ana$constants$SOURCE_DIR,"R","setup_paths.R"))
source(file.path(campa_ana$constants$SOURCE_DIR,"R","mixed_models.R"))
source(file.path(campa_ana$constants$SOURCE_DIR,"R","channels_rename_info.R"))
source(file.path(campa_ana$constants$SOURCE_DIR,"R","io.R"))

# create directory to hold plots
plot_dir <- file.path(campa_ana$constants$SOURCE_DIR,"figures","dot_plots")
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
}

# read metadata files
channels_metadata <- read_csv(file.path(DATA_DIR,"channels_metadata.csv")) %>% 
  select(-1)
wells_metadata <- read_csv(file.path(DATA_DIR,"wells_metadata.csv")) %>% 
  select(-1) %>%
  mutate(treatment=if_else(perturbation_duration %in% c("normal","DMSO-120","DMSO-720"),
                           "Control",
                           perturbation_duration),
         treatment=str_replace(treatment,"-120", " (2.5h)"),
         treatment=str_replace(treatment,"-30", " (1h)"),
         treatment=str_replace(treatment,"-720", " (12.5h)"))

# select required wells and their paths
selected_wells <- wells_metadata %>%
  filter(cell_type=="184A1") %>%
  filter(EU==30) %>%
  filter(!secondary_only) %>%
  select(treatment,well_name) %>%
  left_join(campa_res_dirs,by="well_name") %>%
  left_join(data_dirs,by="well_name")

# load campa results
campa_res <- selected_wells %>%
  mutate(model_res = map(campa_res_dir,~read_csv(file.path(.,"features.csv"),show_col_types = FALSE) %>% select(-1))) %>%
  select(-well_name) %>%
  unnest(model_res)

# load additional data (EU)
cell_features <- read_csv(file=file.path(campa_ana$constants$SOURCE_DIR,"additional_data","184A1_cell_features.csv"))

# keep whole nucleus mean intensities (in "all" CSL)
whole_nucleus_intensities <- campa_res %>%
  filter(cluster=="all") %>% 
  inner_join(cell_features, by = c("mapobject_id", "well_name")) %>%
  rename(`00_EU`=Nuclei_Intensity_mean_00_EU) %>%
  select(-Intensity_sum_22_SE) %>%
  pivot_longer(matches("^\\d{2}_"),names_to="channel",values_to="intensity") 

# find small number of nuclei with zero EU intensity (these are omitted)
whole_nucleus_intensities %>%
  filter(intensity <= 0) %>% 
  nrow()

# Calculate fold changes ----

whole_nucleus_intensity_fold_changes <- map_dfr(
  unique(whole_nucleus_intensities$channel),
  ~fit_mixed_model(
    dat = filter(whole_nucleus_intensities,intensity > 0),
    var = intensity,
    channel_name = .,
    transform = "log",
    random_effect = well_name,
    contrast_var = treatment,
    contrast_var_reference = "Control")
)

# adjust channel names for plotting
whole_nucleus_intensity_fold_changes <- whole_nucleus_intensity_fold_changes %>%
  left_join(rename_channels_for_plotting_dict,
            by=c("channel"="old")) %>%
  select(-channel) %>%
  rename(channel=new)

# Make Heatmap ----

# cluster rows
col_clustering <- whole_nucleus_intensity_fold_changes %>%
  select(treatment,fold_change,channel) %>%
  pivot_wider(names_from = channel,values_from=fold_change) %>%
  column_to_rownames("treatment") %>%
  dist() %>%
  hclust(method = "ward.D2")

# cluster rows
row_clustering <- whole_nucleus_intensity_fold_changes %>%
  select(treatment,fold_change,channel) %>%
  pivot_wider(names_from = channel,values_from=fold_change) %>%
  column_to_rownames("treatment") %>%
  t() %>%
  dist() %>%
  hclust(method = "ward.D2")

col_order <- order.dendrogram(as.dendrogram(col_clustering))
row_order <- order.dendrogram(as.dendrogram(row_clustering))

# define a function to make the bubble plot
make_bubble_plot <- function(dat,row_clustering=NULL,col_clustering=NULL,color_limits=c(-2,2)) {
  if(!is.null(col_clustering)) {
    dat <- dat %>%
      mutate(treatment=factor(treatment,levels=col_clustering))
  }
  if (!is.null(row_clustering)) {
    dat <- dat %>%
      mutate(channel=factor(channel,levels=rev(row_clustering)))
  }
  bubble_plot <- dat %>%
    mutate(p.value.adj.BY = p.adjust(p.value, "BY"),
           log2_fold_change = log2(fold_change)) %>%
    mutate(signif_code = case_when(p.value.adj.BY > 0.05 ~ 0,
                                   p.value.adj.BY > 0.01 ~ 1,
                                   TRUE ~ 2)) %>%
    ggplot(aes(x=treatment,
               y=channel,
               col=log2_fold_change,
               size=signif_code)) + 
    geom_point() +
    ylab("Channel name") +
    xlab("Perturbation") +
    scale_radius(name = expression(italic(p)),
                 breaks = c(0,1,2),
                 labels = c("< 0.05","0.05","< 0.01"),
                 range = c(0.8,2),
                 limits = c(0,2)) + 
    scale_colour_gradient2(name= expression(log[2](`fold-change`)),
                           breaks=c(color_limits[1],0,color_limits[2]),
                           low = "blue",mid="white",high = "red",midpoint = 0,
                           limits=color_limits,oob=scales::squish) +
    theme_bw(base_size = 7) + 
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle=90,hjust=1,vjust=0.5),
          axis.text.y = element_text(angle=0),
          legend.direction = "horizontal",
          legend.title = element_text(margin=margin(0,3,0,3)),
          legend.key.size = unit(0.15,"cm"))
  bubble_plot
}

# make the bubble plot
all_pert_bubble <- make_bubble_plot(whole_nucleus_intensity_fold_changes,
                                    col_clustering = col_clustering$labels[col_order],
                                    row_clustering = row_clustering$labels[row_order])

# Save PDF and PNG output
all_pert_bubble + theme(legend.position = "none")
ggsave(filename = file.path(plot_dir,"whole_nucleus_all_perturbations_comp_control_all_cells_dot_plot_vertical.pdf"),
       width=3.8,height=10,units="cm")
ggsave_cairo(filename = file.path(plot_dir,"whole_nucleus_all_perturbations_comp_control_all_cells_dot_plot_vertical.png"),
             width=3.8,height=10,units="cm",dpi=600)

# Save legend separately
ggpubr::as_ggplot(ggpubr::get_legend(all_pert_bubble, position = NULL))
ggsave(filename = file.path(plot_dir,"whole_nucleus_all_perturbations_comp_control_all_cells_dot_plot_vertical_legend.pdf"),
       width=3.8,height=3,units="cm")

# generate dendrogram as a separate plots
tr <- ggtree(as.dendrogram(col_clustering),right=F,size=0.2)
ggsave(filename = file.path(plot_dir,"whole_nucleus_all_perturbations_comp_control_all_cells_dot_plot_tree_vertical_channels.pdf"),width=1,height=10,units="cm")
tr <- ggtree(as.dendrogram(row_clustering),right=F,size=0.2)
ggsave(filename = file.path(plot_dir,"whole_nucleus_all_perturbations_comp_control_all_cells_dot_plot_tree_vertical_treatments.pdf"),width=1,height=3,units="cm")

# Horizontal 
all_pert_bubble + coord_flip() + theme(legend.position = "none")

# Save PDF and PNG output
ggsave(filename = file.path(plot_dir,"whole_nucleus_all_perturbations_comp_control_all_cells_dot_plot_horizontal.pdf"),
       height=3.8,width=10,units="cm")
ggsave_cairo(filename = file.path(plot_dir,"whole_nucleus_all_perturbations_comp_control_all_cells_dot_plot_horizontal.png"),
             height=3.8,width=10,units="cm",dpi=600)

# generate dendrogram as a separate plots
tr <- ggtree(as.dendrogram(col_clustering),right=T,size=0.2)
ggsave(filename = file.path(plot_dir,"whole_nucleus_all_perturbations_comp_control_all_cells_dot_plot_tree_horizontal_channels.pdf"),width=1,height=10,units="cm")
tr <- ggtree(as.dendrogram(row_clustering),right=T,size=0.2)
ggsave(filename = file.path(plot_dir,"whole_nucleus_all_perturbations_comp_control_all_cells_dot_plot_tree_horizontal_treatments.pdf"),width=1,height=3,units="cm")

# Save legend separately
ggpubr::as_ggplot(ggpubr::get_legend(all_pert_bubble + theme(legend.direction = "vertical"), position = NULL))
ggsave(filename = file.path(plot_dir,"whole_nucleus_all_perturbations_comp_control_all_cells_dot_plot_horizontal_legend.pdf"),
       width=3,height=3.8,units="cm")

# DMSO vs control -----

# compute fold changes
whole_nucleus_intensity_fold_changes_dmso_vs_control <- map_dfr(
  unique(whole_nucleus_intensities$channel),
  ~fit_mixed_model(
    dat = filter(whole_nucleus_intensities,intensity > 0) %>% 
      filter(treatment=="Control") %>%
      select(-treatment) %>%
      mutate(treatment=if_else(perturbation_duration=="normal","Untreated","DMSO")),
    var = intensity,
    channel_name = .,
    transform = "log",
    random_effect = well_name,
    contrast_var = treatment,
    contrast_var_reference = "Untreated")
)

# adjust channel names for plotting
whole_nucleus_intensity_fold_changes_dmso_vs_control <- whole_nucleus_intensity_fold_changes_dmso_vs_control %>%
  left_join(rename_channels_for_plotting_dict,
            by=c("channel"="old")) %>%
  select(-channel) %>%
  rename(channel=new)

dmso_bubble <- make_bubble_plot(whole_nucleus_intensity_fold_changes_dmso_vs_control,
                                row_clustering = row_clustering$labels[row_order])

# Save PDF and PNG output
dmso_bubble + theme(legend.position = "none")
ggsave(filename = file.path(plot_dir,"whole_nucleus_DMSO_comp_untreated_all_cells_dot_plot.pdf"),
       width=3.2,height=10,units="cm")
ggsave_cairo(filename = file.path(plot_dir,"whole_nucleus_DMSO_comp_untreated_all_cells_dot_plot.png"),
             width=3.2,height=10,units="cm",dpi=600)

# Save legend separately
ggpubr::as_ggplot(ggpubr::get_legend(dmso_bubble, position = NULL))
ggsave(filename = file.path(plot_dir,"whole_nucleus_DMSO_comp_untreated_all_cells_dot_plot_vertical_legend.pdf"),
       width=3.8,height=3,units="cm")


