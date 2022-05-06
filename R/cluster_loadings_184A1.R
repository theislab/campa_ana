# clear the workspace
rm(list=ls())

# load required packages
library(tidyverse)
library(patchwork)
library(ComplexHeatmap)

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
source(file.path(campa_ana$constants$SOURCE_DIR,"R","channels_rename_info.R"))
source(file.path(campa_ana$constants$SOURCE_DIR,"R","io.R"))
source(file.path(campa_ana$constants$SOURCE_DIR,"R","colours.R"))

# create directory to hold plots
plot_dir <- file.path(campa_ana$constants$SOURCE_DIR,"figures","cluster_loadings")
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
}

# read cluster annotation
annotation <- read_csv(file.path(campa_ana$constants$SOURCE_DIR,"R","annotation_VAE_all.csv"))
# make a color look-up table
getcolor <- distinct(annotation,cluster_annotation,color) %>% pull(color)
names(getcolor) <- distinct(annotation,cluster_annotation,color) %>% pull(cluster_annotation)

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

# load campa results (per cell)
campa_leiden_res <- selected_wells %>%
  mutate(model_res = map(campa_res_dir,~read_csv(file.path(.,"features.csv"),show_col_types = FALSE) %>% select(-1))) %>%
  select(-well_name) %>%
  unnest(model_res)

# load campa results (per cell) - annotated
campa_csl_res <- selected_wells %>%
  mutate(model_res = map(campa_res_dir,~read_csv(file.path(.,"features_annotation.csv"),show_col_types = FALSE) %>% select(-1))) %>%
  select(-well_name) %>%
  unnest(model_res)

# Make plots ----
# PLOT - show mapping of leiden clusters to compartments using loading matrix
unperturbed_by_leiden_id <- campa_leiden_res %>%
  filter(treatment=="Control") %>%
  # convery to sum intensity to average across all pixels rather than averaging per-cell values
  mutate(across(matches("^\\d{2}_"),~.*size)) %>%
  group_by(cluster) %>%
  summarise(across(matches("^\\d{2}_|size"),sum)) %>%
  mutate(across(matches("^\\d{2}_"),~./size)) %>%
  select(-size) %>%
  filter(cluster!="all") %>%
  column_to_rownames("cluster") %>%
  data.matrix() %>%
  scale() 

order2 <- seriation::seriate(dist(t(unperturbed_by_leiden_id)),method="OLO_ward")

ha = HeatmapAnnotation(foo = anno_block(gp = gpar(
  fill = c(annotation$color[4],
           annotation$color[9],
           annotation$color[10],
           annotation$color[1],
           annotation$color[7],
           annotation$color[5],
           annotation$color[2]))),
  which="row")

column_labels <- structure(rename_channels_for_plotting_dict$new, names = rename_channels_for_plotting_dict$old)

unperturbed_hm_leiden <- ComplexHeatmap::Heatmap(
  unperturbed_by_leiden_id,
  name="Mean\nIntensity\n(z-score)",
  col = cols_pm(3),
  column_labels = column_labels[colnames(unperturbed_by_leiden_id)],
  cluster_columns = as.dendrogram(order2[[1]]),
  row_dend_side = "left",
  row_names_side = "right",
  row_split = annotation$cluster_annotation,
  left_annotation = ha,
  border = "black",
  row_title_gp = gpar(fontsize=8),
  row_title_rot = 0,
  row_names_gp = gpar(fontsize=8),
  column_names_gp = gpar(fontsize=8),
  heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "bold"),
                              labels_gp = gpar(fontsize = 8))
)
draw(unperturbed_hm_leiden)
pdf(file = file.path(plot_dir,"leiden_cluster_loadings_control.pdf"),width=15/cm(1),height=7/cm(1))
unperturbed_hm_leiden
dev.off()
png(file = file.path(plot_dir,"leiden_cluster_loadings_control.png"),width=15/cm(1),height=7/cm(1),units="in",res=600)
unperturbed_hm_leiden
dev.off()

unperturbed_by_csl <- campa_csl_res %>%
  filter(treatment=="Control") %>%
  # convery to sum intensity to average across all pixels rather than averaging per-cell values
  mutate(across(matches("^\\d{2}_"),~.*size)) %>%
  group_by(cluster) %>%
  summarise(across(matches("^\\d{2}_|size"),sum)) %>%
  mutate(across(matches("^\\d{2}_"),~./size)) %>%
  select(-size) %>%
  filter(cluster!="all") %>%
  column_to_rownames("cluster") %>%
  data.matrix() %>%
  scale() 

order1 <- seriation::seriate(dist(unperturbed_by_csl),method="OLO_ward")
order2 <- seriation::seriate(dist(t(unperturbed_by_csl)),method="OLO_ward")

row_colors <- annotation$color
names(row_colors) <- annotation$cluster_annotation

ha = HeatmapAnnotation(foo = rownames(unperturbed_by_csl),
                       col = list(foo = row_colors),
                       which="row",
                       border = T,
                       show_legend = F,
                       show_annotation_name = F)

unperturbed_hm_annotation <- ComplexHeatmap::Heatmap(
  unperturbed_by_csl,
  name="Mean\nIntensity\n(z-score)",
  col = cols_pm(3),
  column_labels = column_labels[colnames(unperturbed_by_csl)],
  cluster_columns = as.dendrogram(order2[[1]]),
  cluster_rows = as.dendrogram(order1[[1]]),
  row_dend_side = "left",
  row_names_side = "right",
  right_annotation = ha,
  border = "black",
  row_title_gp = gpar(fontsize=8),
  row_title_rot = 0,
  row_names_gp = gpar(fontsize=8),
  column_names_gp = gpar(fontsize=8),
  heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "bold"),
                              labels_gp = gpar(fontsize = 8))
)
draw(unperturbed_hm_annotation)
pdf(file = file.path(plot_dir,"csl_cluster_loadings_control.pdf"),width=14/cm(1),height=5.5/cm(1))
draw(unperturbed_hm_annotation)
dev.off()
png(file = file.path(plot_dir,"csl_cluster_loadings_control.png"),width=14/cm(1),height=5.5/cm(1),units="in",res=600)
draw(unperturbed_hm_annotation)
dev.off()

