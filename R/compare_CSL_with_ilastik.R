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
  left_join(data_dirs,by="well_name") %>%
  left_join(ilastik_data_dirs,by="well_name")

# load experimental data
pixels <- selected_wells %>%
  head(1) %>%
  mutate(
    # load names of channels and ilastik channels
    channels = map(data_dir,load_channels),
    ilastik_channels = map(ilastik_data_dir,load_channels),
    # load pixel-level data
    input_data = map2(data_dir,channels,~load_input_data(dirname=.x,channels=.y)),
    ilastik_data = map2(ilastik_data_dir,ilastik_channels,~load_ilastik_data(dirname=.x,channels=.y)),
    # load model clustering output     
    leiden_cluster_id = map(campa_res_dir,~load_pixel_clustering(.,"clustering_res0.5.npy"))) 

ilastik_data_pixels <- pixels %>%
  select(-input_data,-channels,-ilastik_channels,-leiden_cluster_id) %>%
  unnest(ilastik_data) 

campa_pixels <- pixels %>%
  select(-input_data,-channels,-ilastik_channels,-ilastik_data) %>%
  unnest(leiden_cluster_id) %>%
  select(mapobject_id,x,y,cluster_id) %>%
  left_join(annotation,by="cluster_id") %>% 
  mutate(cluster_annotation=factor(cluster_annotation))

comparison_pixels <- left_join(
  ilastik_data_pixels,
  campa_pixels,
  by=c("mapobject_id","x","y")) 

# plot campa clustering
plot_campa <- comparison_pixels %>% 
  filter(mapobject_id==271296) %>%
  ggplot(aes(x=x,y=-y,fill=cluster_annotation)) + 
  geom_raster() +
  scale_fill_manual(values=getcolor[levels(comparison_pixels$cluster_annotation)]) +
  coord_fixed() +
  theme_void(base_size = 7) +
  theme(legend.title=element_blank(),
        legend.direction="vertical",
        legend.position = "bottom",
        legend.key.size = unit(3,"mm"))
plot_campa
ggsave(filename = file.path("PLOTS/CHECK_SEGMENTATION/CAMPA_example_cell.pdf"),
       width=3,height=6,units="cm")
ggsave_cairo(filename = file.path("PLOTS/CHECK_SEGMENTATION/CAMPA_example_cell.png"),
             width=3,height=6,units="cm",dpi=600)

# plot intensity data
plot_intensity_speckles <- pixels %>%
  unnest(input_data) %>%
  filter(mapobject_id==271296) %>%
  pivot_longer(cols=matches('09_SRRM2|15_SON'),names_to="channel",values_to="intensity") %>%
  group_by(channel) %>%
  mutate(intensity = intensity - 105,
         intensity = intensity / quantile(intensity,0.999),
         channel = str_remove(channel,"\\d{2}_")) %>%
  ggplot(aes(x=x,y=-y,fill=intensity)) + 
  facet_wrap(~channel,ncol=1) +
  geom_raster() +
  scale_fill_gradient(low = 'black', high = 'white',limits=c(0,1),oob=scales::oob_squish,breaks=c(0,0.5,1)) +
  coord_fixed() +
  theme_void() 
plot_intensity_speckles

plot_intensity_pml <- pixels %>%
  unnest(input_data) %>%
  filter(mapobject_id==271296) %>%
  pivot_longer(cols=matches('11_PML|20_SP100'),names_to="channel",values_to="intensity") %>%
  group_by(channel) %>%
  mutate(intensity = intensity - 105,
         intensity = intensity / quantile(intensity,0.995),
         channel = str_remove(channel,"\\d{2}_")) %>%
  ggplot(aes(x=x,y=-y,fill=intensity)) + 
  facet_wrap(~channel,ncol=1) +
  geom_raster() +
  scale_fill_gradient(low = 'black', high = 'white',limits=c(0,1),oob=scales::oob_squish,breaks=c(0,0.5,1)) +
  coord_fixed() +
  theme_void() 
plot_intensity_pml

# plot ilastik probability maps 
plot_pm_speckles <- comparison_pixels %>%
  filter(mapobject_id==271296) %>%
  pivot_longer(cols=matches('^\\d{2}'),names_to="ilastik_channel",values_to="probability") %>%
  mutate(ilastik_channel = recode(ilastik_channel,
                                  `09_SRRM2_ILASTIK`="SRRM2-probability",
                                  `11_PML_ILASTIK`="PML-probability",
                                  `15_SON_ILASTIK`="SON-probability",
                                  `20_SP100_ILASTIK`="SP100-probability",)) %>%
  filter(ilastik_channel %in% c("SRRM2-probability","SON-probability")) %>%
  mutate(probability=probability/(2^16-1)) %>%
  ggplot(aes(x=x,y=-y,fill=probability)) + 
  facet_wrap(~ilastik_channel,ncol=1) +
  geom_raster() +
  scale_fill_viridis_c(breaks=c(0,0.5,1)) +
  coord_fixed() +
  theme_void() 
plot_pm_speckles

plot_pm_pml <- comparison_pixels %>%
  filter(mapobject_id==271296) %>%
  pivot_longer(cols=matches('^\\d{2}'),names_to="ilastik_channel",values_to="probability") %>%
  mutate(ilastik_channel = recode(ilastik_channel,
                                  `09_SRRM2_ILASTIK`="SRRM2-probability",
                                  `11_PML_ILASTIK`="PML-probability",
                                  `15_SON_ILASTIK`="SON-probability",
                                  `20_SP100_ILASTIK`="SP100-probability",)) %>%
  filter(ilastik_channel %in% c("PML-probability","SP100-probability")) %>%
  mutate(probability=probability/(2^16-1)) %>%
  ggplot(aes(x=x,y=-y,fill=probability)) + 
  facet_wrap(~ilastik_channel,ncol=1) +
  geom_raster() +
  scale_fill_viridis_c(breaks=c(0,0.5,1)) +
  coord_fixed() +
  theme_void() 
plot_pm_pml

# plot thresholded probability maps
plot_classification_speckles <- comparison_pixels %>%
  filter(mapobject_id==271296) %>%
  mutate(`SRRM2-classifier` = if_else(`09_SRRM2_ILASTIK` > (2^16)*0.95,T,F),
         `SON-classifier` = if_else(`15_SON_ILASTIK` > (2^16)*0.95,T,F),
         `PML-classifier` = if_else(`11_PML_ILASTIK` > (2^16)*0.6,T,F),
         `SP100-classifier` = if_else(`20_SP100_ILASTIK` > (2^16)*0.6,T,F) ) %>%
  pivot_longer(cols=matches('classifier$'),names_to="classifier") %>%
  filter(classifier %in% c("SRRM2-classifier","SON-classifier")) %>%
  ggplot(aes(x=x,y=-y,fill=value)) + 
  facet_wrap(~classifier,ncol=1) +
  geom_raster() +
  scale_fill_grey() +
  coord_fixed() +
  theme_void() +
  theme(legend.title=element_blank())
plot_classification_speckles

plot_classification_pml <- comparison_pixels %>%
  filter(mapobject_id==271296) %>%
  mutate(`SRRM2-classifier` = if_else(`09_SRRM2_ILASTIK` > (2^16)*0.95,T,F),
         `SON-classifier` = if_else(`15_SON_ILASTIK` > (2^16)*0.95,T,F),
         `PML-classifier` = if_else(`11_PML_ILASTIK` > (2^16)*0.6,T,F),
         `SP100-classifier` = if_else(`20_SP100_ILASTIK` > (2^16)*0.6,T,F) ) %>%
  pivot_longer(cols=matches('classifier$'),names_to="classifier") %>%
  filter(classifier %in% c("PML-classifier","SP100-classifier")) %>%
  ggplot(aes(x=x,y=-y,fill=value)) + 
  facet_wrap(~classifier,ncol=1) +
  geom_raster() +
  scale_fill_grey() +
  coord_fixed() +
  theme_void() +
  theme(legend.title=element_blank())
plot_classification_pml

plot_intensity_speckles + 
  plot_pm_speckles + 
  plot_classification_speckles & 
  theme_void(base_size = 7) +
  theme(legend.position = "bottom",
        legend.title=element_blank(),
        legend.key.size = unit(3,"mm"))
ggsave(filename = file.path("PLOTS/CHECK_SEGMENTATION/Nuclear_speckles_intensity_classifiers.pdf"),
       width=8,height=6,units="cm")
ggsave_cairo(filename = file.path("PLOTS/CHECK_SEGMENTATION/Nuclear_speckles_intensity_classifiers.png"),
             width=8,height=6,units="cm",dpi=600)

plot_intensity_pml + 
  plot_pm_pml + 
  plot_classification_pml & 
  theme_void(base_size = 7) +
  theme(legend.position = "bottom",
        legend.title=element_blank(),
        legend.key.size = unit(3,"mm"))
ggsave(filename = file.path("PLOTS/CHECK_SEGMENTATION/PML_bodies_intensity_classifiers.pdf"),
       width=8,height=6,units="cm")
ggsave_cairo(filename = file.path("PLOTS/CHECK_SEGMENTATION/PML_bodies_intensity_classifiers.png"),
             width=8,height=6,units="cm",dpi=600)


# compare clustering for this example cell ----
compare_all_example <- comparison_pixels %>%
  filter(mapobject_id==271296) %>%
  mutate(`SRRM2-classifier` = if_else(`09_SRRM2_ILASTIK` > (2^16)*0.95,T,F),
         `SON-classifier` = if_else(`15_SON_ILASTIK` > (2^16)*0.95,T,F),
         `PML-classifier` = if_else(`11_PML_ILASTIK` > (2^16)*0.6,T,F),
         `SP100-classifier` = if_else(`20_SP100_ILASTIK` > (2^16)*0.6,T,F) ) %>%
  mutate(
    comparison_campa_srrm2 = case_when(
      !`SRRM2-classifier` & cluster_annotation != "Nuclear speckles" ~ "Both false",
      `SRRM2-classifier` & cluster_annotation != "Nuclear speckles" ~ "SRRM2-classifier only",
      !`SRRM2-classifier` & cluster_annotation == "Nuclear speckles" ~ "CAMPA only",
      `SRRM2-classifier` & cluster_annotation == "Nuclear speckles" ~ "Both true"),
    comparison_campa_son = case_when(
      !`SON-classifier` & cluster_annotation != "Nuclear speckles" ~ "Both false",
      `SON-classifier` & cluster_annotation != "Nuclear speckles" ~ "SON-classifier only",
      !`SON-classifier` & cluster_annotation == "Nuclear speckles" ~ "CAMPA only",
      `SON-classifier` & cluster_annotation == "Nuclear speckles" ~ "Both true"),
    comparison_son_srrm2 = case_when(
      !`SON-classifier` & !`SRRM2-classifier` ~ "Both false",
      `SON-classifier` & !`SRRM2-classifier` ~ "SON-classifier only",
      !`SON-classifier` & `SRRM2-classifier` ~ "SRRM2-classifier only",
      `SON-classifier` & `SRRM2-classifier` ~ "Both true"),
    comparison_campa_pml = case_when(
      !`PML-classifier` & cluster_annotation != "PML bodies" ~ "Both false",
      `PML-classifier` & cluster_annotation != "PML bodies" ~ "PML-classifier only",
      !`PML-classifier` & cluster_annotation == "PML bodies" ~ "CAMPA only",
      `PML-classifier` & cluster_annotation == "PML bodies" ~ "Both true"),
    comparison_campa_sp100 = case_when(
      !`SP100-classifier` & cluster_annotation != "PML bodies" ~ "Both false",
      `SP100-classifier` & cluster_annotation != "PML bodies" ~ "SP100-classifier only",
      !`SP100-classifier` & cluster_annotation == "PML bodies" ~ "CAMPA only",
      `SP100-classifier` & cluster_annotation == "PML bodies" ~ "Both true"),
    comparison_sp100_pml = case_when(
      !`SP100-classifier` & !`PML-classifier` ~ "Both false",
      `SP100-classifier` & !`PML-classifier` ~ "SP100-classifier only",
      !`SP100-classifier` & `PML-classifier` ~ "PML-classifier only",
      `SP100-classifier` & `PML-classifier` ~ "Both true")) 

plot_comparison_speckles_example_1 <- compare_all_example %>%
  ggplot(aes(x=x,y=-y,fill=comparison_campa_srrm2)) + 
  geom_raster() +
  scale_fill_manual(values=c("#E3D0D8",'#7FD8BE','#bc2c1a','#315659')) +
  coord_fixed() +
  theme_void() +
  theme(legend.title=element_blank())

plot_comparison_speckles_example_2 <- compare_all_example %>%
  ggplot(aes(x=x,y=-y,fill=comparison_campa_son)) + 
  geom_raster() +
  scale_fill_manual(values=c("#E3D0D8",'#7FD8BE','#bc2c1a',"#8075FF")) +
  coord_fixed() +
  theme_void() +
  theme(legend.title=element_blank())

plot_comparison_speckles_example_3 <- compare_all_example %>%
  ggplot(aes(x=x,y=-y,fill=comparison_son_srrm2)) + 
  geom_raster() +
  scale_fill_manual(values=c("#E3D0D8",'#7FD8BE',"#8075FF",'#315659')) +
  coord_fixed() +
  theme_void() +
  theme(legend.title=element_blank())

plot_comparison_pml_example_1 <- compare_all_example %>%
  ggplot(aes(x=x,y=-y,fill=comparison_campa_pml)) + 
  geom_raster() +
  scale_fill_manual(values=c("#E3D0D8",'#7FD8BE','#bc2c1a','#315659')) +
  coord_fixed() +
  theme_void() +
  theme(legend.title=element_blank())

plot_comparison_pml_example_2 <- compare_all_example %>%
  ggplot(aes(x=x,y=-y,fill=comparison_campa_sp100)) + 
  geom_raster() +
  scale_fill_manual(values=c("#E3D0D8",'#7FD8BE','#bc2c1a',"#8075FF")) +
  coord_fixed() +
  theme_void() +
  theme(legend.title=element_blank())

plot_comparison_pml_example_3 <- compare_all_example %>%
  ggplot(aes(x=x,y=-y,fill=comparison_sp100_pml)) + 
  geom_raster() +
  scale_fill_manual(values=c("#E3D0D8",'#7FD8BE',"#8075FF",'#315659')) +
  coord_fixed() +
  theme_void() +
  theme(legend.title=element_blank())

plot_comparison_speckles_example_1 + 
  plot_comparison_speckles_example_2 + 
  plot_comparison_speckles_example_3 &
  theme_void(base_size = 7) +
  theme(legend.position = "bottom",
        legend.direction = "vertical",
        legend.key.size = unit(3,"mm"),
        legend.title = element_blank()) 
ggsave(filename = file.path("PLOTS/CHECK_SEGMENTATION/segmentation_comparison_speckles.pdf"),
       width=8,height=5,units="cm")
ggsave_cairo(filename = file.path("PLOTS/CHECK_SEGMENTATION/segmentation_comparison_speckles.png"),
             width=8,height=5,units="cm",dpi=600)

plot_comparison_pml_example_1 + 
  plot_comparison_pml_example_2 + 
  plot_comparison_pml_example_3 & 
  theme_void(base_size = 7) +
  theme(legend.position = "bottom",
        legend.direction = "vertical",
        legend.key.size = unit(3,"mm"),
        legend.title = element_blank()) 
ggsave(filename = file.path("PLOTS/CHECK_SEGMENTATION/segmentation_comparison_pml.pdf"),
       width=8,height=5,units="cm")
ggsave_cairo(filename = file.path("PLOTS/CHECK_SEGMENTATION/segmentation_comparison_pml.png"),
             width=8,height=5,units="cm",dpi=600)


# quantitative comparison
compare_pixels_speckles <- function(dat,threshold=0.95) {
  require(yardstick)
  one_hot <- dat %>%
    mutate(`SRRM2-classifier` = as_factor(if_else(`09_SRRM2_ILASTIK` > (2^16)*threshold,T,F)),
           `SON-classifier` = as_factor(if_else(`15_SON_ILASTIK` > (2^16)*threshold,T,F))) %>%
    select(`SRRM2-classifier`,`SON-classifier`,cluster_annotation) %>%
    mutate(campa_nuclear_speckles = as_factor(if_else(cluster_annotation=="Nuclear speckles",T,F)))
  
  multi_metric <- metric_set(accuracy,bal_accuracy,f_meas,sensitivity,specificity)
  
  srrm2_vs_son <- multi_metric(data = one_hot,
                               truth = `SRRM2-classifier`,
                               estimate = `SON-classifier`) %>%
    mutate(comparison = "SON-classifier vs. SRRM2-classifier")
  
  campa_vs_son <- multi_metric(data = one_hot,
                               truth = `SON-classifier`,
                               estimate = campa_nuclear_speckles) %>%
    mutate(comparison = "CAMPA vs. SON-classifier")
  
  campa_vs_srrm2 <- multi_metric(data = one_hot,
                                 truth = `SRRM2-classifier`,
                                 estimate = campa_nuclear_speckles) %>%
    mutate(comparison = "CAMPA vs. SRRM2-classifier")
  
  return(bind_rows(srrm2_vs_son,campa_vs_son,campa_vs_srrm2) %>% mutate(pm_threshold = threshold))
}

compare_pixels_pml <- function(dat,threshold=0.95) {
  require(yardstick)
  one_hot <- dat %>%
    mutate(`PML-classifier` = as_factor(if_else(`11_PML_ILASTIK` > (2^16)*threshold,T,F)),
           `SP100-classifier` = as_factor(if_else(`20_SP100_ILASTIK` > (2^16)*threshold,T,F) )) %>%
    select(`PML-classifier`,`SP100-classifier`,cluster_annotation) %>%
    mutate(campa_pml_bodies = as_factor(if_else(cluster_annotation=="PML bodies",T,F)))
  
  multi_metric <- metric_set(accuracy,bal_accuracy,f_meas,sensitivity,specificity)
  
  sp100_vs_pml <- multi_metric(data = one_hot,
                               truth = `PML-classifier`,
                               estimate = `SP100-classifier`) %>%
    mutate(comparison = "SP100-classifier vs. PML-classifier")
  
  campa_vs_pml <- multi_metric(data = one_hot,
                               truth = `PML-classifier`,
                               estimate = campa_pml_bodies) %>%
    mutate(comparison = "CAMPA vs. PML-classifier")
  
  campa_vs_sp100 <- multi_metric(data = one_hot,
                                 truth = `SP100-classifier`,
                                 estimate = campa_pml_bodies) %>%
    mutate(comparison = "CAMPA vs. SP100-classifier")
  
  return(bind_rows(sp100_vs_pml,campa_vs_pml,campa_vs_sp100) %>% mutate(pm_threshold = threshold))
}


comparison_pixels_sample <- comparison_pixels %>%
  group_by(well_name,treatment) %>%
  sample_frac(0.1) %>%
  ungroup()

# tmp <- comparison_pixels %>%
#   group_by(well_name,treatment) %>%
#   distinct(mapobject_id) %>%
#   sample_frac(0.05) %>%
#   inner_join(comparison_pixels) %>%
#   group_by(mapobject_id) %>%
#   nest() %>%
#   mutate(comp = map(data,~compare_pixels_speckles(.,threshold = 0.95)))
# 
tmp %>% unnest(comp) %>% filter(comparison=="CAMPA vs. SON-classifier") %>%
  filter(.metric=="bal_accuracy") %>%
  arrange(-.estimate)

comp_results_speckles <- map_dfr(c(seq(0.25,0.95,by=0.05),0.99),
                                 ~compare_pixels_speckles(comparison_pixels_sample,threshold = .))

ggplot(comp_results_speckles,aes(x=pm_threshold,y=.estimate,col=comparison)) + 
  geom_line() + 
  facet_wrap(~.metric) +
  coord_cartesian(ylim = c(0,1)) +
  geom_hline(yintercept=0.5,lty="dashed") +
  geom_vline(xintercept=0.95,lty="dotted") +
  theme_bw()

comp_results_pml <- map_dfr(seq(0.25,0.95,by=0.05),
                            ~compare_pixels_pml(comparison_pixels_sample,threshold = .))

ggplot(comp_results_pml,aes(x=pm_threshold,y=.estimate,col=comparison)) + 
  geom_line() + 
  facet_wrap(~.metric) +
  coord_cartesian(ylim = c(0,1)) +
  geom_hline(yintercept=0.5,lty="dashed") +
  geom_vline(xintercept=0.6,lty="dotted") +
  theme_bw()






