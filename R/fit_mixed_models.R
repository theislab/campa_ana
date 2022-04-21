# clear the workspace
rm(list=ls())

# load required packages
library(tidyverse)
library(nlme)
library(emmeans)

# load parallel libraries for multicore processing
library(parallel)
library(foreach)

# set n_cores_override = NULL to use maximum number of cores
n_cores_override <- 12

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

# create directory to hold results of mixed models
model_dir <- file.path(campa_ana$constants$SOURCE_DIR,"mixed_model_results")
if (!dir.exists(model_dir)) {
  dir.create(model_dir)
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
  mutate(model_res = map(campa_res_dir,~read_csv(file.path(.,"features_annotation.csv"),show_col_types = FALSE) %>% select(-1))) %>%
  select(-well_name) %>%
  unnest(model_res) 

# for the purposes of computing statistical differences,
# we exclude a small number of missing compartments
not_identified <- campa_res %>%
  filter(size==0) %>%
  select(treatment,mapobject_id,cluster)

# summarise these
fraction_not_identified <- campa_res %>%
  mutate(zero_size = size==0) %>%
  group_by(treatment,cluster,zero_size) %>%
  count() %>%
  group_by(treatment,cluster) %>%
  mutate(frac_identified = n/sum(n)) %>%
  filter(zero_size == FALSE)

ggplot(fraction_not_identified,aes(x=cluster,y=frac_identified,fill=treatment)) + 
  geom_col(position = "dodge")

# and find a few structures with intensity of zero (problematic for log transformation)
non_positive_intensity <- campa_res %>%
  filter(size!=0) %>%
  pivot_longer(matches("^\\d{2}_"),names_to="channel",values_to="intensity") %>%
  filter(intensity <=0) %>%
  select(treatment,mapobject_id,cluster,channel)

# summarise these
fraction_non_positive_intensity <- campa_res %>%
  filter(size!=0) %>%
  pivot_longer(matches("^\\d{2}_"),names_to="channel",values_to="intensity") %>%
  mutate(non_positive_intensity = intensity <=0) %>%
  group_by(treatment,cluster,non_positive_intensity) %>%
  count() %>%
  group_by(treatment,cluster) %>%
  mutate(frac_positive_intensity = n/sum(n)) %>%
  filter(non_positive_intensity == FALSE)

ggplot(fraction_non_positive_intensity,aes(x=cluster,y=frac_positive_intensity,fill=treatment)) + 
  geom_col(position = "dodge")

# remove
campa_res_exclude_missing <- campa_res %>%
  filter(size>0)

# setup parallel processing
n_cores_available <- detectCores()
n_cores <- ifelse(is.null(n_cores_override),n_cores_available,n_cores_override)
message(paste0('Detected ',n_cores_available,"cores for parallel computation. Using ", n_cores,"."))
cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)

# compute intensity fold-changes across all treatments and channels separately ----

# convert to long format and remove non-positive
campa_res_long_intensities <- campa_res_exclude_missing %>%
  pivot_longer(matches("^\\d{2}_"),names_to="channel",values_to="intensity") %>%
  filter(intensity > 0 ) 

all_treatments <- setdiff(unique(campa_res_long_intensities$treatment),"Control")
all_channels <- unique(campa_res_long_intensities$channel)

system.time(
  intensity_fold_changes_list <- foreach (i=1:length(all_channels)) %dopar% {
    res <- purrr::map_dfr(
      all_treatments,
      ~fit_mixed_model_per_CSL(
        dat = filter(campa_res_long_intensities,treatment %in% c("Control",.)),
        var = intensity,
        channel_name = all_channels[i],
        transform = "log",
        object_id = mapobject_id,
        random_effect = well_name,
        contrast_var = treatment,
        contrast_var_reference = "Control",
        group_var = cluster,
        group_var_reference = "all",
        unnormalised_only = F)
    )
  }
)

# save intensity fold-changes
intensity_fold_changes <- bind_rows(intensity_fold_changes_list)
write_csv(x = intensity_fold_changes,
          file = file.path(model_dir,"184A1_intensity_fold_changes.csv"))

# compute size fold-changes across all treatments separately ----

# extract appropriate data
campa_res_sizes <- campa_res_exclude_missing %>%
  select(-matches("^\\d{2}_"))

# process treatments in parallel
system.time(
  size_fold_changes_list <- foreach (i=1:length(all_treatments)) %dopar%
    fit_mixed_model_per_CSL(
      dat = filter(campa_res_sizes,treatment %in% c("Control",all_treatments[i])),
      var = size,
      transform = "log",
      object_id = mapobject_id,
      random_effect = well_name,
      contrast_var = treatment,
      contrast_var_reference = "Control",
      group_var = cluster,
      group_var_reference = "all",
      unnormalised_only = F)
)

# save intensity fold-changes
size_fold_changes <- bind_rows(size_fold_changes_list)
write_csv(x = size_fold_changes,
          file = file.path(model_dir,"184A1_size_fold_changes.csv"))
