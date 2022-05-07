# clear the workspace
rm(list=ls())

# load required packages
library(tidyverse)
library(patchwork)

# setup script-specific parameters
cell_type <- "HeLa"
experiment_name <- "VAE_SBF2/CondVAE_siRNA-CC"
example_mapobject_ids <- c(200942,256477)

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

# read cluster annotation
annotation <- read_csv(file.path(campa_ana$constants$SOURCE_DIR,"R","annotation_VAE_SBF2.csv"))
# make a color look-up table
getcolor <- distinct(annotation,cluster_annotation,color) %>% pull(color)
names(getcolor) <- distinct(annotation,cluster_annotation,color) %>% pull(cluster_annotation)

# read metadata files
channels_metadata <- read_csv(file.path(DATA_DIR,"channels_metadata.csv")) %>% 
  select(-1)
wells_metadata <- read_csv(file.path(DATA_DIR,"wells_metadata.csv")) %>% 
  select(-1) 

# select required wells and their paths
selected_wells <- wells_metadata %>%
  filter(cell_type=="HeLa") %>%
  filter(EU==30) %>%
  filter(!secondary_only) %>%
  select(siRNA,well_name) %>%
  inner_join(campa_res_dirs,by="well_name") %>%
  left_join(data_dirs,by="well_name")

# find the wells containing the selected mapobject_ids
selected_wells <- selected_wells %>%
  filter(!is.na(data_dir)) %>%
  ungroup() %>%
  mutate(mapobject_id_data = map(data_dir,~load_obj_ids(dirname=.))) %>%
  unnest(mapobject_id_data) %>%
  filter(mapobject_id %in% example_mapobject_ids) %>%
  distinct(siRNA,well_name,data_dir,campa_res_dir) 

# load experimental data
cell_aggregate_data <- selected_wells %>%
  mutate(model_res = map(campa_res_dir,~read_csv(file.path(.,"features_seed3_annotation.csv"),show_col_types = FALSE) %>% select(-1))) %>%
  select(-well_name,-siRNA) %>%
  unnest(model_res) %>%
  mutate(siRNA = factor(siRNA,levels=c("scrambled","SBF2")))


pixels <- selected_wells %>%
  mutate(
    # load names of channels
    channels = map(data_dir,load_channels),
    # load pixel-level data
    input_data = map2(data_dir,channels,~load_input_data(dirname=.x,channels=.y)),
    # load model clustering output     
    leiden_cluster_id = map(campa_res_dir,~load_pixel_clustering(.,"clustering_res0.9_sub-0.33_seed3.npy"))) 

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
  left_join(annotation,by=c("cluster_id"="cluster")) %>% 
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

# calculate size
selected_pixels_outlines %>% 
  filter(mapobject_id==256477) %>%
  ungroup() %>%
  summarise(min_x = min(x), max_x = max(x),
            min_y = min(y), max_y = max(y)) %>%
  mutate(x_size = (max_x-min_x) *6.5/60,
         y_size = (max_y-min_y) *6.5/60)

segmentation_images <- selected_pixels_outlines %>%
  group_by(mapobject_id) %>%
  mutate(x=x-min(x),y=y-min(y)) %>%
  
  ungroup() %>% {
    ggplot(data=.,aes(x=x,y=y,fill=cluster_annotation)) + 
      geom_tile() + 
      geom_tile(data=filter(.,outline==T),fill=grey(0.4)) +
      facet_wrap(~mapobject_id) + 
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
ggsave(segmentation_images,file=file.path(plot_dir,"scrambled_SBF2_examples_full_segmentation.pdf"),width=12,height=7,units="cm")
ggsave_cairo(segmentation_images,file=file.path(plot_dir,"scrambled_SBF2_examples_full_segmentation.png"),width=12,height=7,units="cm",dpi=600)

compartment <- "Nucleolus"
selected_pixels_outlines %>%
  mutate(x=x-min(x),y=y-min(y)) %>%
  ungroup() %>% {
    ggplot(data=.,aes(x=x,y=y,fill=cluster_annotation)) + 
      geom_tile(data=filter(.,cluster_annotation==compartment)) + 
      geom_tile(data=filter(.,outline==T),fill=grey(0.4)) +
      facet_wrap(~mapobject_id) + 
      scale_x_continuous(expand = c(0.05,0.05)) +
      scale_y_continuous(expand = c(0.05,0.05)) +
      theme_void(base_size = 7) +
      coord_fixed() +
      theme(legend.position = "none",
            panel.spacing = unit(0,"mm"),
            strip.text = element_blank()) +
      scale_fill_manual(values = annotation$color[match(c(compartment),annotation$cluster_annotation)])
  }


individual_compartments <- distinct(selected_pixels_outlines,cluster_annotation) %>%
  pull(cluster_annotation)
segmentation_images <- map(
  as.character(individual_compartments),
  function(compartment) {
    selected_pixels_outlines %>%
      mutate(x=x-min(x),y=y-min(y)) %>%
      ungroup() %>% {
        ggplot(data=.,aes(x=x,y=y,fill=cluster_annotation)) + 
          geom_tile(data=filter(.,cluster_annotation==compartment)) + 
          geom_tile(data=filter(.,outline==T),fill=grey(0.4)) +
          facet_wrap(~mapobject_id) + 
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
map2(str_replace(as.character(individual_compartments),"/","_"),
     segmentation_images,
     ~ggsave(.y,file=file.path(plot_dir,paste0("scrambled_SBF2_examples_",.x,"_segmentation.pdf")),width=12,height=7,units="cm"))
map2(str_replace(as.character(individual_compartments),"/","_"),
     segmentation_images,
     ~ggsave_cairo(.y,file=file.path(plot_dir,paste0("scrambled_SBF2_examples_",.x,"_segmentation.png")),width=12,height=7,units="cm",dpi=600))


plot_single_channel <- function(dat,channel_name,percentile_rescaling=0.995) {
  channel_name <- enquo(channel_name)
  channel_image <- dat %>%
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
        geom_tile(data=filter(.,outline==T),fill="#c49a6c",alpha=0.5) +
        facet_wrap(~mapobject_id) + 
        scale_x_continuous(expand = c(.05,.05)) +
        scale_y_continuous(expand = c(.05,.05)) +
        facet_grid(~new) + 
        scale_fill_gradient(low="white",high="black",limits=c(0,1.2),breaks=c(0,0.5,1),oob=scales::squish) +
        theme_void(base_size = 7) +
        coord_fixed() + 
        theme(legend.position = "none")
    }
  ggsave(channel_image,file=file.path(plot_dir,paste0("scrambled_SBF2_examples_",rlang::as_name(channel_name),".pdf")),width=8,height=3,units="cm")
  ggsave_cairo(channel_image,file=file.path(plot_dir,paste0("scrambled_SBF2_examples_",rlang::as_name(channel_name),".png")),width=8,height=3,units="cm",dpi=600)
  
  channel_image
}

plot_single_channel(selected_pixels_outlines,`06_CTNNB1`)
plot_single_channel(selected_pixels_outlines,`08_PXN`)
plot_single_channel(selected_pixels_outlines,`16_GOLGA2`)
plot_single_channel(selected_pixels_outlines,`14_HSPD1`)
plot_single_channel(selected_pixels_outlines,`19_TUBA1A`)

# Plot with outlines of individual bodies
dummy_pixels <- selected_pixels_outlines %>%
  group_by(mapobject_id) %>%
  select(x,y) %>%
  pivot_longer(cols=c(x,y)) %>%
  group_by(mapobject_id,name) %>%
  summarise(min = min(value),
            max = max(value)) %>%
  pivot_wider(values_from = c(min,max))
dummy_pixels <- bind_rows(select(dummy_pixels,mapobject_id,x=min_x,y=min_y),
                          select(dummy_pixels,mapobject_id,x=max_x,y=max_y))


nucleolus_outlines <- selected_pixels_outlines %>%
  filter(cluster_annotation %in% c("Nucleolus")) %>%
  group_by(mapobject_id) %>%
  nest() %>%
  mutate(data = map(data,~erode_clusters(d=.,aggregate_by = cluster_annotation,width=3,remove_small = T))) %>%
  unnest(cols = c(data)) %>%
  select(mapobject_id,x,y) %>%
  bind_rows(dummy_pixels) %>%
  group_by(mapobject_id) %>%
  nest() %>%
  mutate(sparse = map(data,to_sparse)) %>% 
  mutate(sparse = map(sparse,~dilate_sparse(.,width=7))) %>% 
  mutate(outlines = map(sparse,~to_outlines_sparse(.,size=1))) %>%
  select(-data,-sparse) %>%
  unnest(cols=outlines) %>%
  mutate(outline_nucleolus=T,y_new = x - 1, x_new = y - 1) %>%
  select(outline_nucleolus,y=y_new,x=x_new)

cajal_outlines <- selected_pixels_outlines %>%
  filter(cluster_annotation %in% c("Cajal bodies")) %>%
  group_by(mapobject_id) %>%
  nest() %>%
  mutate(data = map(data,~erode_clusters(d=.,aggregate_by = cluster_annotation,width=3,remove_small = T))) %>%
  unnest(cols = c(data)) %>%
  select(mapobject_id,x,y) %>%
  bind_rows(dummy_pixels) %>%
  group_by(mapobject_id) %>%
  nest() %>%
  mutate(sparse = map(data,to_sparse)) %>% 
  mutate(sparse = map(sparse,~dilate_sparse(.,width=7))) %>% 
  mutate(outlines = map(sparse,~to_outlines_sparse(.,size=1))) %>%
  select(-data,-sparse) %>%
  unnest(cols=outlines) %>%
  mutate(outline_cajal=T,y_new = x - 1, x_new = y - 1) %>%
  select(outline_cajal,y=y_new,x=x_new)

p_outlines <- selected_pixels_outlines %>%
  filter(cluster_annotation %in% c("P-bodies")) %>%
  group_by(mapobject_id) %>%
  nest() %>%
  mutate(data = map(data,~erode_clusters(d=.,aggregate_by = cluster_annotation,width=3,remove_small = T))) %>%
  unnest(cols = c(data)) %>%
  select(mapobject_id,x,y) %>%
  bind_rows(dummy_pixels) %>%
  group_by(mapobject_id) %>%
  nest() %>%
  mutate(sparse = map(data,to_sparse)) %>% 
  mutate(sparse = map(sparse,~dilate_sparse(.,width=7))) %>% 
  mutate(outlines = map(sparse,~to_outlines_sparse(.,size=1))) %>%
  select(-data,-sparse) %>%
  unnest(cols=outlines) %>%
  mutate(outline_p=T,y_new = x - 1, x_new = y - 1) %>%
  select(outline_p,y=y_new,x=x_new)

# plot intensity images
NCL<- selected_pixels_outlines %>%
  left_join(nucleolus_outlines) %>%
  ungroup() %>%
  select(-`00_EU`) %>%
  mutate_at(vars(matches("\\d{2}_")),~./ quantile(.,0.999)) %>% {
    ggplot(data=.,aes(x=x,y=y,fill=`21_NCL`)) + 
      geom_tile() + 
      geom_tile(data=filter(.,outline==T),fill="#c49a6c",alpha=0.5) + 
      geom_tile(data=filter(.,outline_nucleolus==T),fill=getcolor["Nucleolus"][[1]]) } +
  coord_fixed() + 
  facet_wrap(~siRNA,ncol=1,strip.position = "right",as.table = F) +
  scale_fill_gradient(low="white",high="black",limits=c(0,1.2),breaks=c(0,0.5,1),oob=scales::squish) +
  theme_void(base_size = 7) +
  theme(legend.position = "none")
NCL

ggsave(plot = NCL,
       filename = file.path(plot_dir,"NCL_outlines.pdf"),
       width=4,height=5.5,units="cm")
ggsave_cairo(plot = NCL,
             filename = file.path(plot_dir,"NCL_outlines.png"),
             width=4,height=5.5,units="cm",dpi=600)

COIL <- selected_pixels_outlines %>%
  left_join(cajal_outlines) %>%
  ungroup() %>%
  select(-`00_EU`) %>%
  mutate_at(vars(matches("\\d{2}_")),~./ quantile(.,0.999)) %>% {
    ggplot(data=.,aes(x=x,y=y,fill=`21_COIL`)) + 
      geom_tile() + 
      geom_tile(data=filter(.,outline==T),fill="#c49a6c",alpha=0.5) + 
      geom_tile(data=filter(.,outline_cajal==T),fill=getcolor["Cajal bodies"][[1]]) } +
  coord_fixed() + 
  facet_wrap(~siRNA,ncol=1,strip.position = "right",as.table = F) +
  scale_fill_gradient(low="white",high="black",limits=c(0,1.2),breaks=c(0,0.5,1),oob=scales::squish) +
  theme_void(base_size = 7) +
  theme(legend.position = "none")
COIL

ggsave(plot = COIL,
       filename = file.path(plot_dir,"COIL_outlines.pdf"),
       width=4,height=5.5,units="cm")
ggsave_cairo(plot = COIL,
             filename = file.path(plot_dir,"COIL_outlines.png"),
             width=4,height=5.5,units="cm",dpi = 600)

DDX6 <- selected_pixels_outlines %>%
  left_join(p_outlines) %>%
  ungroup() %>%
  select(-`00_EU`) %>%
  mutate_at(vars(matches("\\d{2}_")),~./ quantile(.,0.999)) %>% {
    ggplot(data=.,aes(x=x,y=y,fill=`22_DDX6`)) + 
      geom_tile() + 
      geom_tile(data=filter(.,outline==T),fill="#c49a6c",alpha=0.5) + 
      geom_tile(data=filter(.,outline_p==T),fill=getcolor["P-bodies"][[1]]) } +
  coord_fixed() + 
  facet_wrap(~siRNA,ncol=1,strip.position = "right",as.table = F) +
  scale_fill_gradient(low="white",high="black",limits=c(0,1.2),breaks=c(0,0.5,1),oob=scales::squish) +
  theme_void(base_size = 7) +
  theme(legend.position = "none")
DDX6

ggsave(plot = DDX6,
       filename = file.path(plot_dir,"DDX6_outlines.pdf"),
       width=4,height=5.5,units="cm")
ggsave_cairo(plot = DDX6,
             filename = file.path(plot_dir,"DDX6_outlines.png"),
             width=4,height=5.5,units="cm",dpi = 600)


