# clear the workspace
rm(list=ls())

# load required packages
library(tidyverse)
library(patchwork)

# setup script-specific parameters
cell_type <- "184A1"
experiment_name <- "VAE_all/CondVAE_pert-CC"
example_mapobject_id <- 271296

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

# create directory to hold plots
plot_dir <- file.path(campa_ana$constants$SOURCE_DIR,"figures","compare_ilastik_segmentation")
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
  left_join(data_dirs,by="well_name") %>%
  left_join(ilastik_data_dirs,by="well_name")

# load experimental data
set.seed(11)
pixels <- selected_wells %>%
  filter(!is.na(data_dir)) %>%
  group_by(treatment) %>%
  # one well per treatment (otherwise too memory-intensive)
  sample_n(1) %>%
  ungroup() %>%
  mutate(
    # load names of channels and ilastik channels
    channels = map(data_dir,load_channels),
    ilastik_channels = map(ilastik_data_dir,load_channels),
    # load pixel-level data
    input_data = map2(data_dir,channels,~load_input_data(dirname=.x,channels=.y)),
    ilastik_data = map2(ilastik_data_dir,ilastik_channels,~load_ilastik_data(dirname=.x,channels=.y)),
    # load model clustering output     
    leiden_cluster_id = map(campa_res_dir,~load_pixel_clustering(.,"clustering_res0.5.npy"))) 

# get ilastik data only
ilastik_data_pixels <- pixels %>%
  select(-input_data,-channels,-ilastik_channels,-leiden_cluster_id) %>%
  unnest(ilastik_data) %>%
  mutate(mapobject_id = as.integer(mapobject_id)) %>%
  select(-campa_res_dir,-data_dir,-ilastik_data_dir) %>%
  data.table::as.data.table()

# get pixel clustering only
campa_pixels <- pixels %>%
  select(-input_data,-channels,-ilastik_channels,-ilastik_data) %>%
  unnest(leiden_cluster_id) %>%
  select(mapobject_id,x,y,cluster_id) %>%
  left_join(annotation,by="cluster_id") %>% 
  mutate(cluster_annotation=factor(cluster_annotation)) %>%
  mutate(mapobject_id = as.integer(mapobject_id),
         cluster_id = as.integer(cluster_id)) %>%
  data.table::as.data.table() 

# Join is time-consuming in dplyr. Use data.table
comparison_pixels <- data.table::merge.data.table(
  ilastik_data_pixels,
  campa_pixels,
  by=c("mapobject_id","x","y"))

comparison_pixels <- tibble(comparison_pixels)

# plot campa clustering for example cell
plot_campa <- comparison_pixels %>% 
  filter(mapobject_id==example_mapobject_id) %>%
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
ggsave(filename = file.path(plot_dir,"CAMPA_example_cell.pdf"),
       width=3,height=6,units="cm")
ggsave_cairo(filename = file.path(plot_dir,"CAMPA_example_cell.png"),
             width=3,height=6,units="cm",dpi=600)

# plot intensity data for ilastik speckles channels
plot_intensity_speckles <- pixels %>%
  unnest(input_data) %>%
  filter(mapobject_id==example_mapobject_id) %>%
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

# plot intensity data for ilastik pml bodies channels
plot_intensity_pml <- pixels %>%
  unnest(input_data) %>%
  filter(mapobject_id==example_mapobject_id) %>%
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

# plot ilastik probability maps (speckles)
plot_pm_speckles <- comparison_pixels %>%
  filter(mapobject_id==example_mapobject_id) %>%
  pivot_longer(cols=matches('^\\d{2}'),names_to="ilastik_channel",values_to="probability") %>%
  mutate(ilastik_channel = recode(ilastik_channel,
                                  `09_SRRM2_ILASTIK`="SRRM2-probability",
                                  `15_SON_ILASTIK`="SON-probability")) %>%
  filter(ilastik_channel %in% c("SRRM2-probability","SON-probability")) %>%
  mutate(probability=probability/(2^16-1)) %>%
  ggplot(aes(x=x,y=-y,fill=probability)) + 
  facet_wrap(~ilastik_channel,ncol=1) +
  geom_raster() +
  scale_fill_viridis_c(breaks=c(0,0.5,1)) +
  coord_fixed() +
  theme_void() 
plot_pm_speckles

# plot ilastik probability maps (pml bodies)
plot_pm_pml <- comparison_pixels %>%
  filter(mapobject_id==example_mapobject_id) %>%
  pivot_longer(cols=matches('^\\d{2}'),names_to="ilastik_channel",values_to="probability") %>%
  mutate(ilastik_channel = recode(ilastik_channel,
                                  `11_PML_ILASTIK`="PML-probability",
                                  `20_SP100_ILASTIK`="SP100-probability")) %>%
  filter(ilastik_channel %in% c("PML-probability","SP100-probability")) %>%
  mutate(probability=probability/(2^16-1)) %>%
  ggplot(aes(x=x,y=-y,fill=probability)) + 
  facet_wrap(~ilastik_channel,ncol=1) +
  geom_raster() +
  scale_fill_viridis_c(breaks=c(0,0.5,1)) +
  coord_fixed() +
  theme_void() 
plot_pm_pml

# plot thresholded probability maps (speckles)
plot_classification_speckles <- comparison_pixels %>%
  filter(mapobject_id==example_mapobject_id) %>%
  mutate(`SRRM2-classifier` = if_else(`09_SRRM2_ILASTIK` > (2^16)*0.95,T,F),
         `SON-classifier` = if_else(`15_SON_ILASTIK` > (2^16)*0.95,T,F) ) %>%
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

# plot thresholded probability maps (pml)
plot_classification_pml <- comparison_pixels %>%
  filter(mapobject_id==example_mapobject_id) %>%
  mutate(`PML-classifier` = if_else(`11_PML_ILASTIK` > (2^16)*0.6,T,F),
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

# combine plots and save (speckles)
plot_intensity_speckles + 
  plot_pm_speckles + 
  plot_classification_speckles & 
  theme_void(base_size = 7) +
  theme(legend.position = "bottom",
        legend.title=element_blank(),
        legend.key.size = unit(3,"mm"))
ggsave(filename = file.path(plot_dir,"Nuclear_speckles_intensity_classifiers.pdf"),
       width=8,height=6,units="cm")
ggsave_cairo(filename = file.path(plot_dir,"Nuclear_speckles_intensity_classifiers.png"),
             width=8,height=6,units="cm",dpi=600)

# combine plots and save (pml)
plot_intensity_pml + 
  plot_pm_pml + 
  plot_classification_pml & 
  theme_void(base_size = 7) +
  theme(legend.position = "bottom",
        legend.title=element_blank(),
        legend.key.size = unit(3,"mm"))
ggsave(filename = file.path(plot_dir,"PML_bodies_intensity_classifiers.pdf"),
       width=8,height=6,units="cm")
ggsave_cairo(filename = file.path(plot_dir,"PML_bodies_intensity_classifiers.png"),
             width=8,height=6,units="cm",dpi=600)

# compare CAMPA and ilastik clustering for example cell ----
compare_all_example <- comparison_pixels %>%
  filter(mapobject_id==example_mapobject_id) %>%
  mutate(`SRRM2-classifier` = if_else(`09_SRRM2_ILASTIK` > (2^16)*0.95,T,F),
         `SON-classifier` = if_else(`15_SON_ILASTIK` > (2^16)*0.95,T,F),
         `PML-classifier` = if_else(`11_PML_ILASTIK` > (2^16)*0.6,T,F),
         `SP100-classifier` = if_else(`20_SP100_ILASTIK` > (2^16)*0.6,T,F) ) %>%
  mutate(
    comparison_campa_srrm2 = case_when(
      !`SRRM2-classifier` & cluster_annotation != "Nuclear speckles" ~ "Both false",
      `SRRM2-classifier` & cluster_annotation != "Nuclear speckles" ~ "SRRM2-classifier only",
      !`SRRM2-classifier` & cluster_annotation == "Nuclear speckles" ~ "Nucl. speckles CSL only",
      `SRRM2-classifier` & cluster_annotation == "Nuclear speckles" ~ "Both true"),
    comparison_campa_son = case_when(
      !`SON-classifier` & cluster_annotation != "Nuclear speckles" ~ "Both false",
      `SON-classifier` & cluster_annotation != "Nuclear speckles" ~ "SON-classifier only",
      !`SON-classifier` & cluster_annotation == "Nuclear speckles" ~ "Nucl. speckles CSL only",
      `SON-classifier` & cluster_annotation == "Nuclear speckles" ~ "Both true"),
    comparison_son_srrm2 = case_when(
      !`SON-classifier` & !`SRRM2-classifier` ~ "Both false",
      `SON-classifier` & !`SRRM2-classifier` ~ "SON-classifier only",
      !`SON-classifier` & `SRRM2-classifier` ~ "SRRM2-classifier only",
      `SON-classifier` & `SRRM2-classifier` ~ "Both true"),
    comparison_campa_pml = case_when(
      !`PML-classifier` & cluster_annotation != "PML bodies" ~ "Both false",
      `PML-classifier` & cluster_annotation != "PML bodies" ~ "PML-classifier only",
      !`PML-classifier` & cluster_annotation == "PML bodies" ~ "PML-bodies CSL only",
      `PML-classifier` & cluster_annotation == "PML bodies" ~ "Both true"),
    comparison_campa_sp100 = case_when(
      !`SP100-classifier` & cluster_annotation != "PML bodies" ~ "Both false",
      `SP100-classifier` & cluster_annotation != "PML bodies" ~ "SP100-classifier only",
      !`SP100-classifier` & cluster_annotation == "PML bodies" ~ "PML-bodies CSL only",
      `SP100-classifier` & cluster_annotation == "PML bodies" ~ "Both true"),
    comparison_sp100_pml = case_when(
      !`SP100-classifier` & !`PML-classifier` ~ "Both false",
      `SP100-classifier` & !`PML-classifier` ~ "SP100-classifier only",
      !`SP100-classifier` & `PML-classifier` ~ "PML-classifier only",
      `SP100-classifier` & `PML-classifier` ~ "Both true")) 

# compare campa and srrm2
plot_comparison_speckles_example_1 <- compare_all_example %>%
  ggplot(aes(x=x,y=-y,fill=comparison_campa_srrm2)) + 
  geom_raster() +
  scale_fill_manual(values=c("#E3D0D8",'#7FD8BE','#bc2c1a','#315659')) +
  coord_fixed() +
  theme_void() +
  theme(legend.title=element_blank())

# compare campa and son
plot_comparison_speckles_example_2 <- compare_all_example %>%
  ggplot(aes(x=x,y=-y,fill=comparison_campa_son)) + 
  geom_raster() +
  scale_fill_manual(values=c("#E3D0D8",'#7FD8BE','#bc2c1a',"#8075FF")) +
  coord_fixed() +
  theme_void() +
  theme(legend.title=element_blank())

# compare son and srrm2
plot_comparison_speckles_example_3 <- compare_all_example %>%
  ggplot(aes(x=x,y=-y,fill=comparison_son_srrm2)) + 
  geom_raster() +
  scale_fill_manual(values=c("#E3D0D8",'#7FD8BE',"#8075FF",'#315659')) +
  coord_fixed() +
  theme_void() +
  theme(legend.title=element_blank())

# compare campa and pml
plot_comparison_pml_example_1 <- compare_all_example %>%
  ggplot(aes(x=x,y=-y,fill=comparison_campa_pml)) + 
  geom_raster() +
  scale_fill_manual(values=c("#E3D0D8",'#7FD8BE','#bc2c1a','#315659')) +
  coord_fixed() +
  theme_void() +
  theme(legend.title=element_blank())

# compare campa and sp100
plot_comparison_pml_example_2 <- compare_all_example %>%
  ggplot(aes(x=x,y=-y,fill=comparison_campa_sp100)) + 
  geom_raster() +
  scale_fill_manual(values=c("#E3D0D8",'#7FD8BE','#bc2c1a',"#8075FF")) +
  coord_fixed() +
  theme_void() +
  theme(legend.title=element_blank())

# compare sp100 and pml
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
        legend.key.size = unit(2,"mm"),
        legend.title = element_blank()) 
ggsave(filename = file.path(plot_dir,"segmentation_comparison_speckles.pdf"),
       width=8.5,height=5,units="cm")
ggsave_cairo(filename = file.path(plot_dir,"segmentation_comparison_speckles.png"),
             width=8.5,height=5,units="cm",dpi=600)

plot_comparison_pml_example_1 + 
  plot_comparison_pml_example_2 + 
  plot_comparison_pml_example_3 & 
  theme_void(base_size = 7) +
  theme(legend.position = "bottom",
        legend.direction = "vertical",
        legend.key.size = unit(2,"mm"),
        legend.title = element_blank()) 
ggsave(filename = file.path(plot_dir,"segmentation_comparison_pml.pdf"),
       width=8.5,height=5,units="cm")
ggsave_cairo(filename = file.path(plot_dir,"segmentation_comparison_pml.png"),
             width=8.5,height=5,units="cm",dpi=600)


# quantitative comparison function (speckles)
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
  
  son_vs_srrm2 <- multi_metric(data = one_hot,
                               truth = `SON-classifier`,
                               estimate = `SRRM2-classifier`) %>%
    mutate(comparison = "SRRM2-classifier vs. SON-classifier")
  
  campa_vs_son <- multi_metric(data = one_hot,
                               truth = `SON-classifier`,
                               estimate = campa_nuclear_speckles) %>%
    mutate(comparison = "Nuclear speckles CSL vs. SON-classifier")
  
  campa_vs_srrm2 <- multi_metric(data = one_hot,
                                 truth = `SRRM2-classifier`,
                                 estimate = campa_nuclear_speckles) %>%
    mutate(comparison = "Nuclear speckles CSL vs. SRRM2-classifier")
  
  return(bind_rows(srrm2_vs_son,son_vs_srrm2,campa_vs_son,campa_vs_srrm2) %>% mutate(pm_threshold = threshold))
}

# quantitative comparison function (pml)
compare_pixels_pml <- function(dat,threshold=0.95) {
  require(yardstick)
  one_hot <- dat %>%
    mutate(`PML-classifier` = as_factor(if_else(`11_PML_ILASTIK` > (2^16)*threshold,T,F)),
           `SP100-classifier` = as_factor(if_else(`20_SP100_ILASTIK` > (2^16)*threshold,T,F) )) %>%
    select(`PML-classifier`,`SP100-classifier`,cluster_annotation) %>%
    mutate(campa_pml_bodies = as_factor(if_else(cluster_annotation=="PML bodies",T,F)))
  
  multi_metric <- metric_set(accuracy,bal_accuracy,f_meas,sensitivity,specificity)
  
  pml_vs_sp100 <- multi_metric(data = one_hot,
                               truth = `PML-classifier`,
                               estimate = `SP100-classifier`) %>%
    mutate(comparison = "SP100-classifier vs. PML-classifier")
  
  sp100_vs_pml <- multi_metric(data = one_hot,
                               truth = `SP100-classifier`,
                               estimate = `PML-classifier`) %>%
    mutate(comparison = "PML-classifier vs. SP100-classifier")
  
  campa_vs_pml <- multi_metric(data = one_hot,
                               truth = `PML-classifier`,
                               estimate = campa_pml_bodies) %>%
    mutate(comparison = "PML-bodies CSL vs. PML-classifier")
  
  campa_vs_sp100 <- multi_metric(data = one_hot,
                                 truth = `SP100-classifier`,
                                 estimate = campa_pml_bodies) %>%
    mutate(comparison = "PML-bodies CSL vs. SP100-classifier")
  
  return(bind_rows(pml_vs_sp100,sp100_vs_pml,campa_vs_pml,campa_vs_sp100) %>% mutate(pm_threshold = threshold))
}

# Sub-sample for quantitative comparison
comparison_pixels_sample <- comparison_pixels %>%
  group_by(well_name,treatment) %>%
  sample_frac(0.05) %>%
  ungroup()

# compute comparison over many different values of the threshold
comp_results_speckles <- map_dfr(c(seq(0.25,0.95,by=0.05),0.99),
                                 ~compare_pixels_speckles(comparison_pixels_sample,threshold = .))

comp_results_speckles %>%
  filter(.metric %in% c("f_meas","sens","spec")) %>%
  mutate(.metric = recode(.metric,"f_meas"="F1-score","sens"="Sensitivity","spec"="Specificity")) %>%
  ggplot(aes(x=pm_threshold,y=.estimate,col=comparison)) + 
  geom_line(size=0.3) + 
  facet_wrap(~.metric) +
  coord_cartesian(ylim = c(0,1)) +
  geom_hline(yintercept=0.5,lty="75",size=0.2) +
  geom_hline(yintercept=1.0,lty="solid",size=0.2) +
  geom_vline(xintercept=0.95,lty="23",size=0.2) +
  theme_bw(base_size = 7) +
  theme(panel.grid = element_blank(),
        legend.key.size = unit(0.3,"cm"),
        legend.title=element_blank(),
        legend.position = "bottom",
        legend.direction = "vertical") +
  scale_colour_brewer(palette = "Dark2") +
  scale_y_continuous(name="Metric value",breaks=c(0,0.5,1)) +
  scale_x_continuous(name="Probability map threshold",breaks=c(0,0.5,1),limits = c(0.4,1))
ggsave(filename = file.path(plot_dir,"CAMPA_compared_to_ilastik_metrics_speckles.pdf"),
       width=7,height=5.5,units="cm")
ggsave_cairo(filename = file.path(plot_dir,"CAMPA_compared_to_ilastik_metrics_speckles.png"),
             width=7,height=5.5,units="cm",dpi=600)

comp_results_pml <- map_dfr(c(seq(0.25,0.95,by=0.05),0.99),
                            ~compare_pixels_pml(comparison_pixels_sample,threshold = .))

comp_results_pml %>%
  filter(.metric %in% c("f_meas","sens","spec")) %>%
  mutate(.metric = recode(.metric,"f_meas"="F1-score","sens"="Sensitivity","spec"="Specificity")) %>%
  ggplot(aes(x=pm_threshold,y=.estimate,col=comparison)) + 
  geom_line(size=0.3) + 
  facet_wrap(~.metric) +
  coord_cartesian(ylim = c(0,1)) +
  geom_hline(yintercept=0.5,lty="75",size=0.2) +
  geom_hline(yintercept=1.0,lty="solid",size=0.2) +
  geom_vline(xintercept=0.6,lty="23",size=0.2) +
  theme_bw(base_size = 7) +
  theme(panel.grid = element_blank(),
        legend.key.size = unit(0.3,"cm"),
        legend.title=element_blank(),
        legend.position = "bottom",
        legend.direction = "vertical") +
  scale_colour_brewer(palette = "Dark2") +
  scale_y_continuous(name="Metric value",breaks=c(0,0.5,1)) +
  scale_x_continuous(name="Probability map threshold",breaks=c(0,0.5,1),limits = c(0.4,1))
ggsave(filename = file.path(plot_dir,"CAMPA_compared_to_ilastik_metrics_PML.pdf"),
       width=7,height=5.5,units="cm")
ggsave_cairo(filename = file.path(plot_dir,"CAMPA_compared_to_ilastik_metrics_PML.png"),
             width=7,height=5.5,units="cm",dpi=600)


# now for the given threshold, calculate across all conditions
all_comparisons_fixed_threshold <- comparison_pixels %>%
  group_by(well_name,treatment) %>%
  nest() %>%
  mutate(speckles = map(data,~compare_pixels_speckles(.,threshold=0.95)),
         pml = map(data,~compare_pixels_pml(.,threshold=0.6)))

# speckles
speckles_all_comparison <- all_comparisons_fixed_threshold %>%
  select(-data) %>%
  unnest(cols=speckles) %>%
  filter(.metric=="f_meas") 

speckles_all_comparison_summary <- speckles_all_comparison %>%
  group_by(comparison) %>%
  summarise(mean=mean(.estimate),
            sd=sd(.estimate))

speckles_all_comparison %>%
  select(treatment,comparison,.estimate) %>%
  ggplot(aes(x=comparison,y=.estimate)) + 
  geom_col(data=speckles_all_comparison_summary,
           aes(y=mean),
           fill="grey70",size=0.25) +
  geom_errorbar(data=speckles_all_comparison_summary,
                aes(ymin=mean-sd,ymax=mean+sd,y=mean),size=0.25,width=0.2) +
  #geom_jitter(pch=21,fill="white",size=0.5,stroke=0.2) + 
  geom_hline(yintercept=0.5,lty="75",size=0.2) +
  geom_hline(yintercept=1.0,lty="solid",size=0.2) +
  ylab("F1-score") +
  theme_bw(base_size = 7) +
  coord_flip() +
  theme(panel.grid.minor = element_blank(),
        axis.title.y=element_blank())
ggsave(filename = file.path(plot_dir,"F1_speckles.pdf"),
       width=7,height=2.5,units="cm")
ggsave_cairo(filename = file.path(plot_dir,"F1_speckles.png"),
             width=7,height=2.5,units="cm",dpi=600)

speckles_all_comparison_summary

speckles_all_comparison %>%
  ungroup() %>%
  select(-well_name) %>%
  select(treatment,comparison,.estimate) %>%
  filter(comparison!="SRRM2-classifier vs. SON-classifier") %>%
  mutate(label = sprintf("%0.3f",.estimate)) %>%
  ggplot(aes(x=treatment,y=comparison,fill=.estimate)) + 
  geom_tile(col="white") + 
  geom_text(aes(label=label),col="white",size=5/.pt) +
  scale_fill_distiller(name="F1-score",palette = "Reds",type = "seq",direction = 1,limits=c(0,1),breaks=c(0,0.5,1)) +
  theme_bw(base_size = 7) +
  theme(axis.title = element_blank(),
        #aspect.ratio = 1,
        axis.text.x = element_text(angle=90,hjust=1,vjust=0.5),
        legend.key.size = unit(2,"mm"))
ggsave(filename = file.path(plot_dir,"F1_speckles_hm.pdf"),
       width=10,height=3.5,units="cm")
ggsave_cairo(filename = file.path(plot_dir,"F1_speckles_hm.png"),
             width=10,height=3.5,units="cm",dpi=600)
 
# pml
pml_all_comparison <- all_comparisons_fixed_threshold %>%
  select(-data) %>%
  unnest(cols=pml) %>%
  filter(.metric=="f_meas") 

pml_all_comparison_summary <- pml_all_comparison %>%
  group_by(comparison) %>%
  summarise(mean=mean(.estimate),
            sd=sd(.estimate))

pml_all_comparison %>%
  select(treatment,comparison,.estimate) %>%
  ggplot(aes(x=comparison,y=.estimate)) + 
  geom_col(data=pml_all_comparison_summary,
           aes(y=mean),
           fill="grey70",size=0.25) +
  geom_errorbar(data=pml_all_comparison_summary,
             aes(ymin=mean-sd,ymax=mean+sd,y=mean),size=0.25,width=0.2) +
  #geom_jitter(pch=21,fill="white",size=0.5,stroke=0.2) + 
  geom_hline(yintercept=0.5,lty="75",size=0.2) +
  geom_hline(yintercept=1.0,lty="solid",size=0.2) +
  ylab("F1-score") +
  theme_bw(base_size = 7) +
  coord_flip() +
  theme(panel.grid.minor = element_blank(),
        axis.title.y=element_blank())
ggsave(filename = file.path(plot_dir,"F1_pml.pdf"),
       width=8,height=2.5,units="cm")
ggsave_cairo(filename = file.path(plot_dir,"F1_pml.png"),
             width=8,height=2.5,units="cm",dpi=600)

pml_all_comparison %>%
  ungroup() %>%
  select(-well_name) %>%
  select(treatment,comparison,.estimate) %>%
  filter(comparison!="SP100-classifier vs. PML-classifier") %>%
  mutate(label = sprintf("%0.3f",.estimate)) %>%
  ggplot(aes(x=treatment,y=comparison,fill=.estimate)) + 
  geom_tile(col="white") + 
  geom_text(aes(label=label),col="white",size=5/.pt) +
  scale_fill_distiller(name="F1-score",palette = "Reds",type = "seq",direction = 1,limits=c(0,1),breaks=c(0,0.5,1)) +
  theme_bw(base_size = 7) +
  theme(axis.title = element_blank(),
        #aspect.ratio = 1,
        axis.text.x = element_text(angle=90,hjust=1,vjust=0.5),
        legend.key.size = unit(2,"mm"))
ggsave(filename = file.path(plot_dir,"F1_pml_hm.pdf"),
       width=9.5,height=3.5,units="cm")
ggsave_cairo(filename = file.path(plot_dir,"F1_pml_hm.png"),
             width=9.5,height=3.5,units="cm",dpi=600)
