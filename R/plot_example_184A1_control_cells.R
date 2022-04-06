# clear the workspace
rm(list=ls())

# load required packages
library(tidyverse)
library(patchwork)

# setup script-specific parameters
cell_type <- "184A1"
experiment_name <- "VAE_all/CondVAE_pert-CC"
example_mapobject_ids <- c(345916,231622,256562)

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
plot_dir <- file.path(campa_ana$constants$SOURCE_DIR,"figures","example_control_cells")
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

# find the wells containing the selected mapobject_ids
selected_wells <- selected_wells %>%
  filter(!is.na(data_dir)) %>%
  ungroup() %>%
  mutate(mapobject_id_data = map(data_dir,~load_obj_ids(dirname=.))) %>%
  unnest(mapobject_id_data) %>%
  filter(mapobject_id %in% example_mapobject_ids) %>%
  distinct(treatment,well_name,data_dir,campa_res_dir) 

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
selected_pixels <- pivot_longer(selected_pixels,cols = matches("\\d{2}_"),values_to = "intensity") %>%
  left_join(select(channels_metadata,name,mean_background),by="name") %>%
  mutate(intensity = intensity - mean_background) %>%
  select(-mean_background) %>%
  pivot_wider(values_from = intensity)
  
#G1 = 345916
#S = 231622,
#G2 = 256562

# calculate size
selected_pixels_outlines %>% 
  filter(mapobject_id==256562) %>%
  ungroup() %>%
  summarise(min_x = min(x), max_x = max(x),
            min_y = min(y), max_y = max(y)) %>%
  mutate(x_size = (max_x-min_x) *6.5/60,
         y_size = (max_y-min_y) *6.5/60)

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
  group_by(mapobject_id) %>%
  mutate(x=x-min(x),y=y-min(y)) %>%
  
  # HACK!
  mutate(x = if_else(cell_cycle=="S",x+20,x)) %>%
  
  ungroup() %>% {
    ggplot(data=.,aes(x=x,y=y,fill=cluster_annotation)) + 
      geom_tile() + 
      geom_tile(data=filter(.,outline==T),fill=grey(0.4)) +
      facet_wrap(~cell_cycle) + 
      scale_x_continuous(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0)) +
      theme_void(base_size = 7) +
      coord_fixed() + 
      theme(legend.position = "none",
            panel.spacing = unit(0,"mm"),
            strip.text = element_blank()) +
      scale_fill_manual(values = annotation$color[match(levels(selected_pixels_outlines$cluster_annotation),annotation$cluster_annotation)])
  }
segmentation_images
ggsave(segmentation_images,file=file.path(plot_dir,"unperturbed_examples_full_segmentation.pdf"),width=10,height=3,units="cm")
ggsave_cairo(segmentation_images,file=file.path(plot_dir,"unperturbed_examples_full_segmentation.png"),width=10,height=3,units="cm",dpi=600)

individual_compartments <- c("Nucleoplasm","Cajal bodies","Nuclear speckles","PML bodies","Nucleolus","Nuclear periphery")
segmentation_images <- map(
  individual_compartments,
  function(compartment) {
    selected_pixels_outlines %>%
      filter(mapobject_id==256562) %>%
      mutate(x=x-min(x),y=y-min(y)) %>%
      ungroup() %>% {
        ggplot(data=.,aes(x=x,y=y,fill=cluster_annotation)) + 
          geom_tile(data=filter(.,cluster_annotation==compartment)) + 
          geom_tile(data=filter(.,outline==T),fill=grey(0.4)) +
          scale_x_continuous(expand = c(0.05,0.05)) +
          scale_y_continuous(expand = c(0.05,0.05)) +
          theme_void(base_size = 7) +
          coord_fixed() +
          theme(legend.position = "none",
                panel.spacing = unit(0,"mm"),
                strip.text = element_blank()) +
          scale_fill_manual(values = annotation$color[match(c(compartment),annotation$cluster_annotation)])
      }
  })
segmentation_images
map2(individual_compartments,segmentation_images,
     ~ggsave(.y,file=file.path(plot_dir,paste0("unperturbed_examples_256562_",.x,"_segmentation.pdf")),width=1.65,height=1.65,units="cm"))
map2(individual_compartments,segmentation_images,
     ~ggsave_cairo(.y,file=file.path(plot_dir,paste0("unperturbed_examples_256562_",.x,"_segmentation.png")),width=1.65,height=1.65,units="cm",dpi=600))


plot_single_channel <- function(dat,channel_name,percentile_rescaling=0.997) {
  channel_name <- enquo(channel_name)
  channel_image <- dat %>%
    filter(mapobject_id==256562) %>%
    select(mapobject_id,y,x,cell_cycle,outline,
           !! channel_name) %>%
    ungroup() %>%
    mutate(x=x-min(x),y=y-min(y)) %>%
    ungroup() %>%
    pivot_longer(matches("\\d{2}_"),names_to="channel") %>%
    left_join(rename_channels_for_plotting_dict,by=c("channel"="old")) %>%
    group_by(channel) %>% 
    mutate(value = value / quantile(value,percentile_rescaling)) %>% {
      ggplot(data=.,aes(x=x,y=y,fill=value)) + 
        geom_tile() +
        geom_tile(data=filter(.,outline==T),fill=grey(0.4)) +
        scale_x_continuous(expand = c(.05,.05)) +
        scale_y_continuous(expand = c(.05,.05)) +
        facet_grid(~new) + 
        scale_fill_gradient(low = "white",high="black",limits=c(0,1.2),oob=scales::oob_squish) +
        #scale_fill_viridis_c(option ="D",limits=c(0,1.2),oob=scales::oob_squish) +
        theme_void(base_size = 7) +
        coord_fixed() + 
        theme(legend.position = "none")
    }
  ggsave(channel_image,file=file.path(plot_dir,paste0("unperturbed_examples_256562_",rlang::as_name(channel_name),".pdf")),width=1.8,height=1.8,units="cm")
  ggsave_cairo(channel_image,file=file.path(plot_dir,paste0("unperturbed_examples_256562_",rlang::as_name(channel_name),".png")),width=1.8,height=1.8,units="cm",dpi=600)
  
  channel_image
}

plot_single_channel(selected_pixels_outlines,`13_POL2RA_pS5`)
plot_single_channel(selected_pixels_outlines,`21_NCL`)
plot_single_channel(selected_pixels_outlines,`18_NONO`)
plot_single_channel(selected_pixels_outlines,`10_H3K27ac`,percentile_rescaling = 0.99)
plot_single_channel(selected_pixels_outlines,`05_Sm`,percentile_rescaling = 0.99)
plot_single_channel(selected_pixels_outlines,`21_COIL`,percentile_rescaling = 0.9995)
plot_single_channel(selected_pixels_outlines,`11_PML`,percentile_rescaling = 0.999)
plot_single_channel(selected_pixels_outlines,`09_SRRM2`,percentile_rescaling = 0.998)
plot_single_channel(selected_pixels_outlines,`20_ALYREF`)
plot_single_channel(selected_pixels_outlines,`19_KPNA1_MAX`,percentile_rescaling = 0.9995)
plot_single_channel(selected_pixels_outlines,`16_H3`)
plot_single_channel(selected_pixels_outlines,`08_H3K4me3`)

# all cell-cycle phases
individual_compartments <- c("Nucleoplasm","Cajal bodies","Nuclear speckles","PML bodies","Nucleolus")
segmentation_images <- map(
  individual_compartments,
  function(compartment) {
    selected_pixels_outlines %>%
      group_by(mapobject_id) %>%
      mutate(x=x-min(x),y=y-min(y)) %>%
      
      # HACK!
      mutate(x = if_else(cell_cycle=="S",x+20,x)) %>%
      
      ungroup() %>% {
        ggplot(data=.,aes(x=x,y=y,fill=cluster_annotation)) + 
          geom_tile(data=filter(.,cluster_annotation==compartment)) + 
          geom_tile(data=filter(.,outline==T),fill=grey(0.4)) +
          scale_x_continuous(expand = c(0,0)) +
          scale_y_continuous(expand = c(0,0)) +
          facet_wrap(~cell_cycle) + 
          theme_void(base_size = 7) +
          coord_fixed() + 
          theme(legend.position = "none",
                panel.spacing = unit(0,"mm"),
                strip.text = element_blank()) +
          scale_fill_manual(values = annotation$color[match(c(compartment),annotation$cluster_annotation)])
      }
  })
segmentation_images
map2(individual_compartments,segmentation_images,
     ~ggsave(.y,file=file.path(plot_dir,paste0("unperturbed_examples_",.x,"_segmentation.pdf")),width=10,height=3,units="cm"))
map2(individual_compartments,segmentation_images,
     ~ggsave_cairo(.y,file=file.path(plot_dir,paste0("unperturbed_examples_",.x,"_segmentation.png")),width=10,height=3,units="cm",dpi=600))

channel_images <- selected_pixels_outlines %>%
  select(mapobject_id,y,x,cell_cycle,outline,`12_RB1_pS807_S811`,`14_PCNA`) %>%
  ungroup() %>%
  group_by(mapobject_id) %>%
  mutate(x=x-min(x),y=y-min(y)) %>%
  
  # HACK!
  mutate(x = if_else(cell_cycle=="S",x+20,x)) %>%
  
  ungroup() %>%
  pivot_longer(matches("\\d{2}_"),names_to="channel") %>%
  left_join(rename_channels_for_plotting_dict,by=c("channel"="old")) %>%
  group_by(channel) %>% 
  mutate(value = value / quantile(value,0.997)) %>% {
    ggplot(data=.,aes(x=x,y=y,fill=value)) + 
      geom_tile() +
      geom_tile(data=filter(.,outline==T),fill=grey(0.4)) +
      scale_x_continuous(expand = c(.05,.05)) +
      scale_y_continuous(expand = c(.05,.05)) +
      facet_grid(new~cell_cycle) + 
      scale_fill_gradient(low = "white",high="black",limits=c(0,1.2),oob=scales::oob_squish) +
      #scale_fill_viridis_c(option ="D",limits=c(0,1.2),oob=scales::oob_squish) +
      theme_void(base_size = 7) +
      coord_fixed() +
      theme(legend.position = "none")
  }
channel_images
ggsave(channel_images,file=file.path(plot_dir,"unperturbed_examples_channel_images.pdf"),width=6,height=5,units="cm")
ggsave_cairo(channel_images,file=file.path(plot_dir,"unperturbed_examples_channel_images.png"),width=6,height=5,units="cm",dpi=600)
