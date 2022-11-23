# clear the workspace
rm(list=ls())

# load required packages
library(tidyverse)
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
source(file.path(campa_ana$constants$SOURCE_DIR,"R","channels_rename_info.R"))
source(file.path(campa_ana$constants$SOURCE_DIR,"R","io.R"))
source(file.path(campa_ana$constants$SOURCE_DIR,"R","morphology.R"))

# create directory to hold plots
plot_dir <- file.path(campa_ana$constants$SOURCE_DIR,"figures","example_cells")
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
}

current_resolution <- 1.2

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
                           "Unperturbed",
                           perturbation_duration),
         treatment=str_replace(treatment,"-120", " (2.5h)"),
         treatment=str_replace(treatment,"-30", " (1h)"),
         treatment=str_replace(treatment,"-720", " (12.5h)"))

# select required wells and their paths
selected_wells <- wells_metadata %>%
  filter(cell_type=="184A1" & perturbation_duration %in% c("normal","TSA-30")) %>%
  filter(EU==30) %>%
  filter(!secondary_only) %>%
  select(treatment,well_name) %>%
  left_join(campa_res_dirs,by="well_name") %>%
  left_join(data_dirs,by="well_name") 

# find the wells containing the selected mapobject_ids
set.seed(3)
example_mapobject_ids <- selected_wells %>%
  filter(!is.na(data_dir)) %>%
  group_by(treatment) %>%
  top_n(1) %>%
  ungroup() %>%
  mutate(mapobject_id_data = map(data_dir,~load_obj_ids(dirname=.))) %>%
  unnest(mapobject_id_data) %>%
  distinct(treatment,well_name,data_dir,campa_res_dir,mapobject_id) %>%
  group_by(treatment,well_name,data_dir,campa_res_dir) %>%
  sample_n(10) %>%
  ungroup()

selected_wells <- selected_wells %>%
  inner_join(distinct(example_mapobject_ids,well_name)) %>%
  distinct(treatment,well_name,data_dir,campa_res_dir) 

example_mapobject_ids <- pull(example_mapobject_ids,mapobject_id)

# load experimental data
cell_aggregate_data <- selected_wells %>%
  mutate(model_res = map(campa_res_dir,~read_csv(file.path(.,"features_annotation.csv"),show_col_types = FALSE) %>% select(-1))) %>%
  select(-well_name) %>%
  unnest(model_res)

pixels <- selected_wells %>%
  mutate(
    # load names of channels
    channels = map(data_dir,load_channels),
    # load pixel-level data
    input_data = map2(data_dir,channels,~load_input_data(dirname=.x,channels=.y)),
    # load model clustering output     
    leiden_cluster_id = map(campa_res_dir,~load_pixel_clustering(.,"clustering_res0.5.npy")))  

# get input data only
input_data_pixels <- pixels %>%
  select(-channels,-leiden_cluster_id) %>%
  unnest(input_data) %>%
  filter(mapobject_id %in% example_mapobject_ids) %>%
  mutate(mapobject_id = as.integer(mapobject_id)) %>%
  select(-campa_res_dir,-data_dir) %>%
  data.table::as.data.table()

# get pixel clustering only
campa_pixels <- pixels %>%
  select(-input_data,-channels) %>%
  unnest(leiden_cluster_id) %>%
  filter(mapobject_id %in% example_mapobject_ids) %>%
  select(mapobject_id,x,y,cluster_id) %>%
  left_join(annotation,by="cluster_id") %>% 
  mutate(cluster_annotation=factor(cluster_annotation)) %>%
  mutate(mapobject_id = as.integer(mapobject_id),
         cluster_id = as.integer(cluster_id)) %>%
  data.table::as.data.table() 

# Join is time-consuming in dplyr. Use data.table
selected_pixels <- data.table::merge.data.table(
  input_data_pixels,
  campa_pixels,
  by=c("mapobject_id","x","y")) %>%
  tibble() %>%
  left_join(select(cell_aggregate_data,-matches("\\d{2}_"))) %>%
  mutate(cluster_annotation = factor(cluster_annotation),
         cell_cycle = factor(cell_cycle,levels=c("G1","S","G2"))) 

# subtract background
selected_pixels <- selected_pixels %>%
  select(-matches("\\d{2}_(?!H3K27)",perl = T)) %>%
  pivot_longer(cols = matches("\\d{2}_"),values_to = "intensity") %>%
  left_join(select(channels_metadata,name,mean_background),by="name") %>%
  mutate(intensity = intensity - mean_background) %>%
  select(-mean_background) %>%
  pivot_wider(values_from = intensity)

selected_pixels_outlines <- selected_pixels %>%
  group_by(mapobject_id) %>%
  nest() %>%
  mutate(sparse = map(data,to_sparse)) %>%
  mutate(outlines = map(sparse,~to_outlines_sparse(.,size=2))) %>%
  select(-data,-sparse) %>%
  unnest(cols=outlines) %>%
  mutate(outline=T,y_new = x - 1, x_new = y - 1) %>%
  select(outline,y=y_new,x=x_new) %>%
  right_join(selected_pixels)

segmentation_images <- selected_pixels_outlines %>%
  filter(mapobject_id %in% distinct(selected_pixels_outlines,mapobject_id)$mapobject_id[1:8]) %>%
  group_by(mapobject_id) %>%
  mutate(x=x-min(x),y=y-min(y)) %>%
  ungroup() %>% {
    ggplot(data=.,aes(x=x,y=y,fill=cluster_annotation)) + 
      geom_tile() + 
      geom_tile(data=filter(.,outline==T),fill=grey(0.4)) +
      facet_wrap(treatment~mapobject_id,nrow=1) + 
      scale_x_continuous(expand = c(.05,.05)) +
      scale_y_continuous(expand = c(.05,.05)) +
      theme_void(base_size = 7) +
      coord_fixed() + 
      theme(legend.position = "none",
           panel.spacing = unit(0,"mm"),
           strip.text = element_blank())  +
      scale_fill_manual(values = annotation$color[match(levels(as.factor(selected_pixels_outlines$cluster_annotation)),annotation$cluster_annotation)])
  }
segmentation_images
ggsave(segmentation_images,file=file.path(plot_dir,"tsa_versus_control_examples_full_segmentation.pdf"),width=16,height=3.5,units="cm")
ggsave_cairo(segmentation_images,file=file.path(plot_dir,"tsa_versus_control_examples_full_segmentation.png"),width=16,height=3.5,units="cm",dpi=600)


