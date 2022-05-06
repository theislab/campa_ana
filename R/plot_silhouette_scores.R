# clear the workspace
rm(list=ls())

# load required packages
library(tidyverse)
library(patchwork)
library(ComplexHeatmap)

# setup script-specific parameters
cell_type <- "184A1"
experiment_name <- "VAE_all/CondVAE_pert-CC"
pixel_size=6.5/60

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
plot_dir <- file.path(campa_ana$constants$SOURCE_DIR,"figures","silhouette")
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
}

treatment_names <- c('AZD4573 (1h)', 'AZD4573 (2.5h)', 'CX5461 (2.5h)', 'Meayamycin (12.5h)',
                     'TSA (1h)', 'Triptolide (2.5h)', 'Unperturbed')
siscore_co_occ <- read_csv(file.path(campa_ana$constants$SOURCE_DIR,"R","siscore_co_occ.csv")) %>% select(-1) %>% 
  bind_cols(treatment_name = treatment_names)
siscore_mpp <- read_csv(file.path(campa_ana$constants$SOURCE_DIR,"R","siscore_mpp.csv")) %>% select(-1) %>% 
  bind_cols(treatment_name = treatment_names)
siscore_intensity <- read_csv(file.path(campa_ana$constants$SOURCE_DIR,"R","siscore_intensity.csv")) %>% select(-1) %>% 
  bind_cols(treatment_name = treatment_names)

names(siscore_co_occ) <- c(treatment_names,"treatment_name")
names(siscore_mpp) <- c(treatment_names,"treatment_name")
names(siscore_intensity) <- c(treatment_names,"treatment_name")

bind_rows(
  list(pivot_longer(siscore_co_occ,-treatment_name,names_to = "treatment_name2",values_to = "siscore") %>%
         mutate(feature_set="Spatial co-occurrence"),
       pivot_longer(siscore_mpp,-treatment_name,names_to = "treatment_name2",values_to = "siscore") %>%
         mutate(feature_set="Nuclear intensities"),
       pivot_longer(siscore_intensity,-treatment_name,names_to = "treatment_name2",values_to = "siscore") %>%
         mutate(feature_set="CSL intensities"))
) %>%
  filter(treatment_name=="Unperturbed") %>%
  filter(treatment_name2 %in% c("TSA (1h)","CX5461 (2.5h)")) %>%
  mutate(treatment_name2 = str_replace_all(treatment_name2," ","\n")) %>%
  mutate(feature_set=factor(feature_set,levels=c("Nuclear intensities","CSL intensities","Spatial co-occurrence"))) %>%
  ggplot(aes(x=treatment_name2,y=siscore,fill=feature_set)) + 
  geom_col(position = position_dodge(0.7),col="black",size=0.1,width=0.7) +
  theme_bw(base_size = 7) +
  geom_hline(size=0.2,yintercept=0) +
  scale_y_continuous(name = "Difference from unperturbed\n(Silhouette score)", breaks=scales::pretty_breaks(n=5)) +
  scale_fill_brewer(palette="Blues") +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.key.size = unit(3,"mm"),
        legend.position = "bottom",
        legend.direction = "vertical")
ggsave(filename = file.path(plot_dir,"silhouette_scores_selected.pdf"),width=4,height=6,units = "cm")
ggsave_cairo(filename = file.path(plot_dir,"silhouette_scores_selected.pdf"),width=4,height=6,units = "cm")

