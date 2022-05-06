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
# create directory to hold plots
plot_dir <- file.path(campa_ana$constants$SOURCE_DIR,"figures","EU_heterogeneity")
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
}

# read metadata files
channels_metadata <- read_csv(file.path(DATA_DIR,"channels_metadata.csv")) %>% 
  select(-1)
wells_metadata <- read_csv(file.path(DATA_DIR,"wells_metadata.csv")) %>% 
  select(-1) 

cell_features <- read_csv(file=file.path(campa_ana$constants$SOURCE_DIR,"additional_data","HeLa_cell_features.csv"))

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
campa_res <- bind_rows(campa_res,campa_res_nucleus_or_cytoplasm) %>% 
  inner_join(cell_features, by = c("mapobject_id", "well_name")) 

# check cell counts match
length(unique(campa_res_nucleus_or_cytoplasm$mapobject_id))
length(unique(campa_res$mapobject_id))
length(unique(cell_features$mapobject_id))

# keep scrambled cells only
campa_res <- filter(campa_res,siRNA=="scrambled")

# Bin cells by transcription rate ----
cutoffs <- campa_res %>%
  filter(cluster=="all") %>%
  ungroup() %>%
  summarise(extreme_lower_EU = quantile(Nuclei_Intensity_mean_00_EU,0.005),
            lower_EU = quantile(Nuclei_Intensity_mean_00_EU,0.25), # lower quartile
            upper_EU = quantile(Nuclei_Intensity_mean_00_EU,0.75), # upper quartile
            extreme_upper_EU = quantile(Nuclei_Intensity_mean_00_EU,0.995))

campa_res %>%
  filter(cluster=="all") %>%
  ggplot(aes(x=Nuclei_Intensity_mean_00_EU,group=well_name)) +
  geom_rect(data=cutoffs, xmin = cutoffs$extreme_lower_EU,xmax=cutoffs$lower_EU,ymin=-Inf,ymax=Inf,col=NA,fill="blue",alpha=0.3,inherit.aes = F) +
  geom_rect(data=cutoffs, xmin = cutoffs$upper_EU,xmax=cutoffs$extreme_upper_EU,ymin=-Inf,ymax=Inf,col=NA,fill="red",alpha=0.3,inherit.aes = F) +
  geom_density(size=0.2) +
  theme_bw(base_size = 7) +
  ylab("Density") +
  xlab("Mean nuclear EU") +
  scale_x_continuous(breaks = c(cutoffs$extreme_lower_EU,cutoffs$lower_EU,cutoffs$upper_EU,cutoffs$extreme_upper_EU)) +
  coord_cartesian(xlim=c(0,NA),ylim=c(0,0.004),expand=c(0,0)) +
  theme(panel.grid = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(angle = 90,hjust=1,vjust=0.5),
        axis.text.y=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(size=0.2),
        plot.margin = margin(5,10,5,5))
ggsave(filename = file.path(plot_dir,"binning_EU_distribution.pdf"),width=3.5,height=2.5,units="cm")
ggsave_cairo(filename = file.path(plot_dir,"binning_EU_distribution.png"),width=3.5,height=2.5,units="cm",dpi=600)

# apply cutoffs
campa_res <- campa_res %>%
  mutate(EU_bin = case_when(
    Nuclei_Intensity_mean_00_EU > cutoffs$extreme_lower_EU & Nuclei_Intensity_mean_00_EU < cutoffs$lower_EU ~ "lower",
    Nuclei_Intensity_mean_00_EU < cutoffs$extreme_upper_EU & Nuclei_Intensity_mean_00_EU > cutoffs$upper_EU ~ "upper",
    TRUE ~ NA_character_))

# plot binned data
campa_res %>%
  filter(cluster=="all") %>%
  ggplot(aes(x=Nuclei_Intensity_mean_00_EU,col=EU_bin)) +
  geom_density() +
  theme_bw()

# check frequencies of bins by cell cycle stage
campa_res %>%
  filter(cluster=="all") %>%
  group_by(cell_cycle,EU_bin) %>%
  count() %>%
  group_by(EU_bin) %>%
  mutate(freq = n/sum(n))

# exclude non-binned cells
campa_res_bins_only <- campa_res %>%
  filter(!is.na(EU_bin))

# for the purposes of computing statistical differences,
# we exclude a small number of missing compartments
not_identified <- campa_res_bins_only %>%
  filter(size==0) %>%
  select(siRNA,mapobject_id,cluster)

# summarise these
fraction_not_identified <- campa_res_bins_only %>%
  mutate(zero_size = size==0) %>%
  group_by(siRNA,cluster,zero_size) %>%
  count() %>%
  group_by(siRNA,cluster) %>%
  mutate(frac_identified = n/sum(n)) %>%
  filter(zero_size == FALSE)

ggplot(fraction_not_identified,aes(x=cluster,y=frac_identified,fill=siRNA)) + 
  geom_col(position = "dodge")

# and find a few structures with intensity of zero (problematic for log transformation)
non_positive_intensity <- campa_res_bins_only %>%
  filter(size!=0) %>%
  pivot_longer(matches("^\\d{2}_"),names_to="channel",values_to="intensity") %>%
  filter(intensity <=0) %>%
  select(siRNA,mapobject_id,cluster,channel)

# summarise these
fraction_non_positive_intensity <- campa_res_bins_only %>%
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
campa_res_exclude_missing <- campa_res_bins_only %>%
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
               cluster %in% c("P-bodies","all")) %>%
    mutate(EU_bin=factor(EU_bin)),
  var = intensity,
  channel_name = "14_PCNA",
  transform = "log",
  object_id = mapobject_id,
  random_effect = well_name,
  contrast_var = EU_bin,
  contrast_var_reference = "lower",
  group_var = cluster,
  group_var_reference = "all",
  unnormalised_only = T)

# High versus low transcription ----

# parallelise across channels
intensity_fold_changes_list <- foreach (i=1:length(all_channels)) %dopar% {
  # use purrr to map across clusters (individually compared to "all")
  res <- purrr::map_dfr(
    all_clusters,
    ~fit_mixed_model_per_CSL(
      dat = filter(campa_res_long_intensities, cluster %in% c(.,"all")),
      var = intensity,
      channel_name = all_channels[i],
      transform = "log",
      object_id = mapobject_id,
      random_effect = well_name,
      contrast_var = EU_bin,
      contrast_var_reference = "lower",
      group_var = cluster,
      group_var_reference = "all",
      unnormalised_only = F)
  )
}

# save intensity fold-changes
intensity_fold_changes <- dplyr::bind_rows(intensity_fold_changes_list)
write_csv(x = intensity_fold_changes,
          file = file.path(model_dir,"HeLa_intensity_fold_changes_EU_bin.csv"))

# compute size fold-changes across all treatments separately ----

# extract appropriate data
campa_res_sizes <- campa_res_exclude_missing %>%
  filter(cluster!="Antibody aggregates") %>%
  select(-matches("^\\d{2}_")) %>%
  inner_join(subset_objects)

# use purrr to map across clusters (individually compared to "Nucleus (combined)")
size_fold_changes <- purrr::map_dfr(
  all_clusters,
  ~fit_mixed_model_per_CSL(
    dat = filter(campa_res_sizes,cluster %in% c(.,"all")),
    var = size,
    transform = "log",
    object_id = mapobject_id,
    random_effect = well_name,
    contrast_var = EU_bin,
    contrast_var_reference = "lower",
    group_var = cluster,
    group_var_reference = "all",
    unnormalised_only = F)
)

# save size fold-changes
write_csv(x = size_fold_changes,
          file = file.path(model_dir,"HeLa_size_fold_changes_EU_bin.csv"))

stopCluster(cl)



