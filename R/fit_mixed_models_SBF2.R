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
n_cores_override <- 14

# setup script-specific parameters
cell_type <- "HeLa"
experiment_name <- "VAE_SBF2/CondVAE_siRNA-CC"

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
  select(-1) 

# select required wells and their paths
selected_wells <- wells_metadata %>%
  filter(cell_type=="HeLa") %>%
  filter(EU==30) %>%
  filter(!secondary_only) %>%
  select(siRNA,well_name) %>%
  inner_join(campa_res_dirs,by="well_name") %>%
  left_join(data_dirs,by="well_name")

# load campa results
campa_res <- selected_wells %>%
  mutate(model_res = map(campa_res_dir,~read_csv(file.path(.,"features_seed3_annotation.csv"),show_col_types = FALSE) %>% select(-1))) %>%
  select(-well_name,-siRNA) %>%
  unnest(model_res) %>%
  mutate(siRNA = factor(siRNA,levels=c("scrambled","SBF2")))

# generate nuclear and cytoplasmic compartments
nuclear_CSLs <- c(
  "Cajal bodies",             
  "Nuclear periphery",        
  "Nuclear speckles",         
  "Nucleolus",                
  "Nucleoplasm",
  "PML bodies")
cytoplasm_CSLs <- setdiff(unique(campa_res$cluster),c(nuclear_CSLs,"all"))

campa_res_nucleus_or_cytoplasm <- campa_res %>%
  filter(cluster!="all") %>%
  mutate(nucleus_or_cytoplasm = if_else(cluster %in% nuclear_CSLs,"Nucleus (combined)","Cytoplasm (combined)")) %>%
  mutate_at(vars(matches("\\d{2}_")),~.*size) %>%
  group_by(mapobject_id,well_name,siRNA,TR,cell_cycle,nucleus_or_cytoplasm) %>%
  summarise_at(vars(matches("\\d{2}_|size")),~sum(.,na.rm=T)) %>%
  mutate_at(vars(matches("\\d{2}_")),~./size) %>%
  ungroup() %>%
  rename(cluster = nucleus_or_cytoplasm)

# add aggregated nucleus and cytoplasm CSL into the same dataframe
campa_res <- bind_rows(campa_res,campa_res_nucleus_or_cytoplasm) 

# check cell counts match
length(unique(campa_res_nucleus_or_cytoplasm$mapobject_id))
length(unique(campa_res$mapobject_id))

# for the purposes of computing statistical differences,
# we exclude a small number of missing compartments
not_identified <- campa_res %>%
  filter(size==0) %>%
  select(siRNA,mapobject_id,cluster)

# summarise these
fraction_not_identified <- campa_res %>%
  mutate(zero_size = size==0) %>%
  group_by(siRNA,cluster,zero_size) %>%
  count() %>%
  group_by(siRNA,cluster) %>%
  mutate(frac_identified = n/sum(n)) %>%
  filter(zero_size == FALSE)

ggplot(fraction_not_identified,aes(x=cluster,y=frac_identified,fill=siRNA)) + 
  geom_col(position = "dodge")

# and find a few structures with intensity of zero (problematic for log transformation)
non_positive_intensity <- campa_res %>%
  filter(size!=0) %>%
  pivot_longer(matches("^\\d{2}_"),names_to="channel",values_to="intensity") %>%
  filter(intensity <=0) %>%
  select(siRNA,mapobject_id,cluster,channel)

# summarise these
fraction_non_positive_intensity <- campa_res %>%
  filter(size!=0) %>%
  pivot_longer(matches("^\\d{2}_"),names_to="channel",values_to="intensity") %>%
  mutate(non_positive_intensity = intensity <=0) %>%
  group_by(siRNA,cluster,non_positive_intensity) %>%
  count() %>%
  group_by(siRNA,cluster) %>%
  mutate(frac_positive_intensity = n/sum(n)) %>%
  filter(non_positive_intensity == FALSE)

ggplot(fraction_non_positive_intensity,aes(x=cluster,y=frac_positive_intensity,fill=siRNA)) + 
  geom_col(position = "dodge")

# remove
campa_res_exclude_missing <- campa_res %>%
  filter(size>0)

# setup parallel processing
n_cores_available <- detectCores()
n_cores <- ifelse(is.null(n_cores_override),n_cores_available,n_cores_override)
message(paste0('Detected ',n_cores_available," cores for parallel computation. Using ", n_cores,"."))
cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)

# compute intensity fold-changes across all treatments and channels separately ----
subset_objects <- campa_res_exclude_missing %>%
  distinct(mapobject_id) 

# convert to long format and remove non-positive
campa_res_long_intensities <- campa_res_exclude_missing %>%
  filter(cluster!="Antibody aggregates") %>%
  pivot_longer(matches("^\\d{2}_"),names_to="channel",values_to="intensity") %>%
  filter(intensity > 0 ) %>%
  inner_join(subset_objects)

all_channels <- unique(campa_res_long_intensities$channel)
all_clusters <- setdiff(unique(campa_res_long_intensities$cluster),c("all"))

# # example
fit_mixed_model_per_CSL(
  dat = filter(campa_res_long_intensities,
               cluster %in% c("Nuclear periphery","Nucleus (combined)") &
                 siRNA %in% c("scrambled","SBF2")),
  var = intensity,
  channel_name = "03_CDK9",
  transform = "log",
  object_id = mapobject_id,
  random_effect = well_name,
  contrast_var = siRNA,
  contrast_var_reference = "scrambled",
  group_var = cluster,
  group_var_reference = "Nucleus (combined)",
  unnormalised_only = F)

# Fit nuclear CSLs relative to nucleus fold-changes

# parallelise across channels
nuclear_intensity_fold_changes_list <- foreach (i=1:length(all_channels)) %dopar% {
  # use purrr to map across clusters (individually compared to "Nucleus (combined)")
  res <- purrr::map_dfr(
    nuclear_CSLs,
    ~fit_mixed_model_per_CSL(
      dat = filter(campa_res_long_intensities, cluster %in% c(.,"Nucleus (combined)")),
      var = intensity,
      channel_name = all_channels[i],
      transform = "log",
      object_id = mapobject_id,
      random_effect = well_name,
      contrast_var = siRNA,
      contrast_var_reference = "scrambled",
      group_var = cluster,
      group_var_reference = "Nucleus (combined)",
      unnormalised_only = F)
  )
}

# Fit cytoplasmic CSLs relative to cytoplasmic fold-changes

cytoplasm_intensity_fold_changes_list <- foreach (i=1:length(all_channels)) %dopar% {
  # use purrr to map across clusters (individually compared to "Cytoplasm (combined)")
  res <- purrr::map_dfr(
    setdiff(cytoplasm_CSLs,c("Antibody aggregates")),
    ~fit_mixed_model_per_CSL(
      dat = filter(campa_res_long_intensities, cluster %in% c(.,"Cytoplasm (combined)")),
      var = intensity,
      channel_name = all_channels[i],
      transform = "log",
      object_id = mapobject_id,
      random_effect = well_name,
      contrast_var = siRNA,
      contrast_var_reference = "scrambled",
      group_var = cluster,
      group_var_reference = "Cytoplasm (combined)",
      unnormalised_only = F)
  )
}

# save intensity fold-changes
intensity_fold_changes <- dplyr::bind_rows(nuclear_intensity_fold_changes_list,cytoplasm_intensity_fold_changes_list)
write_csv(x = intensity_fold_changes,
          file = file.path(model_dir,"HeLa_intensity_fold_changes_SBF2.csv"))

# compute size fold-changes across all treatments separately ----

# extract appropriate data
campa_res_sizes <- campa_res_exclude_missing %>%
  filter(cluster!="Antibody aggregates") %>%
  select(-matches("^\\d{2}_")) %>%
  inner_join(subset_objects)

# use purrr to map across clusters (individually compared to "Nucleus (combined)")
nuclear_size_fold_changes <- purrr::map_dfr(
  nuclear_CSLs,
  ~fit_mixed_model_per_CSL(
    dat = filter(campa_res_sizes,cluster %in% c(.,"Nucleus (combined)")),
    var = size,
    transform = "log",
    object_id = mapobject_id,
    random_effect = well_name,
    contrast_var = siRNA,
    contrast_var_reference = "scrambled",
    group_var = cluster,
    group_var_reference = "Nucleus (combined)",
    unnormalised_only = F)
)

cytoplasm_size_fold_changes <- purrr::map_dfr(
    setdiff(cytoplasm_CSLs,c("Antibody aggregates")),
  ~fit_mixed_model_per_CSL(
    dat = filter(campa_res_sizes,cluster %in% c(.,"Cytoplasm (combined)")),
    var = size,
    transform = "log",
    object_id = mapobject_id,
    random_effect = well_name,
    contrast_var = siRNA,
    contrast_var_reference = "scrambled",
    group_var = cluster,
    group_var_reference = "Cytoplasm (combined)",
    unnormalised_only = F)
)

# save size fold-changes
size_fold_changes <- dplyr::bind_rows(nuclear_size_fold_changes,cytoplasm_size_fold_changes)
write_csv(x = size_fold_changes,
          file = file.path(model_dir,"HeLa_size_fold_changes_SBF2.csv"))

stopCluster(cl)



