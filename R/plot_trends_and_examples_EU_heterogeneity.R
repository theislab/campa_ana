# clear the workspace
rm(list=ls())

# load required packages
library(tidyverse)

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
source(file.path(campa_ana$constants$SOURCE_DIR,"R","channels_rename_info.R"))
source(file.path(campa_ana$constants$SOURCE_DIR,"R","io.R"))
source(file.path(campa_ana$constants$SOURCE_DIR,"R","mixed_models.R"))
source(file.path(campa_ana$constants$SOURCE_DIR,"R","morphology.R"))

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# read cluster annotation
annotation <- read_csv(file.path(campa_ana$constants$SOURCE_DIR,"R","annotation_VAE_SBF2.csv"))
# make a color look-up table
getcolor <- distinct(annotation,cluster_annotation,color) %>% pull(color)
names(getcolor) <- distinct(annotation,cluster_annotation,color) %>% pull(cluster_annotation)

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


library(boot)
bootstrap_mean_ci <- function(x) {
  out <- boot(data = x,
              R = 500,
              statistic = function(dat, idx) mean(dat[idx], na.rm = TRUE))
  ci <- boot.ci(out, type = "perc",conf = 0.95)
  return(tibble(lower = ci$percent[4],upper = ci$percent[5]))
}

bootstrap_cor_ci <- function(x,y) {
  out <- boot(data = data.frame(x,y),
              R = 500,
              statistic = function(dat, idx) cor(dat[idx,"x"],dat[idx,"y"]))
  ci <- boot.ci(out, type = "perc",conf = 0.95)
  return(tibble(lower = ci$percent[4],upper = ci$percent[5]))
}

bootstrap_cor_all <- function(x,y) {
  out <- boot(data = data.frame(x,y),
              R = 500,
              statistic = function(dat, idx) cor(dat[idx,"x"],dat[idx,"y"]))
  return(tibble(cor = out$t))
}

all_data_EU_trend <- campa_res %>%
  filter(cluster  %in% c("PML bodies")) %>%
  select(matches("PML|SP100|mapobject_id$|Nuclei_Intensity_mean_00_EU")) %>%
  pivot_longer(-c(mapobject_id,Nuclei_Intensity_mean_00_EU),names_to="channel") %>%
  ungroup() %>%
  mutate(EU_bin = cut_interval(Nuclei_Intensity_mean_00_EU,n=30,labels=F))

p_counts <- all_data_EU_trend %>%
  distinct(mapobject_id,EU_bin) %>%
  group_by(EU_bin) %>%
  count() %>%
  filter(n>10) %>%
  ggplot(aes(x=EU_bin-1,y=n)) +
  geom_col(fill="grey70") +
  scale_y_continuous(breaks = c(0,200)) +
  coord_cartesian(ylim=c(0,NA),xlim=c(0,NA),expand = c(0,0)) +
  ylab("Number\nof cells") +
  theme_bw(base_size = 7) +
  theme(panel.border = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.y = element_line(size=0.2),
        axis.line.x = element_line(size=0.2),
        panel.spacing = unit(1,"mm"),
        panel.grid = element_blank())
p_counts

normalisation <- all_data_EU_trend %>%
  group_by(channel) %>%
  summarise(normaliser = mean(value))

p_trends <- all_data_EU_trend %>%
  group_by(channel,EU_bin) %>%
  count() %>%
  filter(n>1) %>%
  inner_join(all_data_EU_trend) %>%
  group_by(channel,EU_bin) %>%
  summarise(n = n(),
            mean = mean(value),
            mean_ci = list(bootstrap_mean_ci(value))) %>%
  filter(n > 10) %>%
  unnest_wider(mean_ci) %>%
  ungroup() %>%
  left_join(normalisation) %>%
  mutate(mean = mean/normaliser,
         lower = lower/normaliser,
         upper = upper/normaliser) %>%
  left_join(rename_channels_for_plotting_dict,by=c("channel"="old")) %>%
  ggplot(aes(x=EU_bin-1,y=mean,col=new)) +
  geom_errorbar(aes(ymin=lower,ymax=upper),width=0.2,size=0.2) +
  geom_smooth(data=all_data_EU_trend %>%
                group_by(channel) %>%
                filter(Nuclei_Intensity_mean_00_EU < quantile(Nuclei_Intensity_mean_00_EU,0.995) &
                         Nuclei_Intensity_mean_00_EU > quantile(Nuclei_Intensity_mean_00_EU,0.005)) %>%
                left_join(rename_channels_for_plotting_dict,by=c("channel"="old")) %>%
                left_join(normalisation),
              aes(y=value/normaliser,fill=new),size=0.2,alpha=0.2) +
  geom_point(pch=21,fill="white",size=0.5) +
  expand_limits(x=0,y=0) +
  scale_color_manual(values=c(cbPalette[7],cbPalette[4])) +
  scale_fill_manual(values=c(cbPalette[7],cbPalette[4])) +
  scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
  scale_x_continuous(limits=c(0,21)) +
  coord_cartesian(ylim=c(0,2),xlim=c(0,NA),expand = c(0,0)) +
  ylab("Mean intensity\n(PML bodies)") +
  xlab("Mean nuclear EU (binned)") +
  theme_bw(base_size = 7) +
  theme(panel.spacing = unit(1,"mm"),
        panel.grid = element_blank(),
        legend.title=element_blank(),
        legend.key.size = unit(2,"mm"))
p_trends

combined_plot <- p_counts / p_trends + plot_layout(heights = c(1,4))
combined_plot
ggsave(plot = combined_plot, filename = file.path(plot_dir,"PML_SP100_trend_EU.pdf"),width=5.5,height=3.8,units="cm")
ggsave_cairo(plot = combined_plot, filename = file.path(plot_dir,"PML_SP100_trend_EU.png"),width=5.5,height=3.8,units="cm",dpi=600)


combined_plot <- p_counts / (p_trends +  
                               scale_color_manual(values=c(cbPalette[3],cbPalette[8])) +
                               scale_fill_manual(values=c(cbPalette[3],cbPalette[8]))) + plot_layout(heights = c(1,4))
combined_plot
ggsave(plot = combined_plot, filename = file.path(plot_dir,"PML_SP100_trend_EU_CM.pdf"),width=5.5,height=3.8,units="cm")
ggsave_cairo(plot = combined_plot, filename = file.path(plot_dir,"PML_SP100_trend_EU_CM.png"),width=5.5,height=3.8,units="cm",dpi=600)

# calculate cell cycle specific correlations with EU
correlations_with_EU <- campa_res %>%
  filter(cluster=="PML bodies") %>%
  select(mapobject_id,`11_PML`,`20_SP100`,Nuclei_Intensity_mean_00_EU,cell_cycle) %>%
  group_by(cell_cycle) %>%
  summarise(cor_PML_EU = cor(`11_PML`,Nuclei_Intensity_mean_00_EU),
            cor_SP100_EU = cor(`20_SP100`,Nuclei_Intensity_mean_00_EU),
            cor_PML_EU_ci = list(bootstrap_cor_ci(`11_PML`,Nuclei_Intensity_mean_00_EU)),
            cor_SP100_EU_ci = list(bootstrap_cor_ci(`20_SP100`,Nuclei_Intensity_mean_00_EU))) %>%
  unnest(col = c(cor_PML_EU_ci,cor_SP100_EU_ci),names_sep="_")

correlations_with_EU_summary <- correlations_with_EU %>%
  pivot_longer(-cell_cycle) %>%
  mutate(channel = if_else(grepl("PML",name),"PML","SP100"),
         statistic = case_when(grepl("lower",name) ~ "lower",
                               grepl("upper",name) ~ "upper",
                               TRUE ~ "cor")) %>%
  select(-name) %>%
  pivot_wider(names_from = statistic) %>%
  mutate(cell_cycle = factor(cell_cycle,levels=c("G1","S","G2"))) 

correlations_with_EU_summary %>%
  ggplot(aes(x=channel,y=cor,fill=cell_cycle)) +
  geom_col(position = position_dodge(0.7),width=0.7,col="black",size=0.2) +
  geom_errorbar(aes(ymin=lower,ymax=upper),position = position_dodge(0.7),size=0.2,width=0.2) +
  geom_hline(yintercept=0,size=0.2) +
  scale_y_continuous("Correlation with EU",limits=c(-.4,.4),breaks=c(-.4,0,.4)) +
  scale_fill_brewer(name="Cell cycle\nstage", palette = "Greys") +
  theme_bw(base_size = 7) +
  #ggtitle("PML bodies") +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_text(size=6),
        #axis.text.x = element_text(angle=90,hjust=1,vjust=0.5),
        legend.key.size = unit(2,"mm"),
        legend.position = "right",
        legend.margin = margin(-4,0,-4,-4),
        legend.direction = "vertical")
ggsave(file.path(plot_dir,"PML_SP100_EU_barplot.pdf"),width=3.5,height=3,units="cm")

# show distributions of bootstraps
correlations_with_EU_all <- campa_res %>%
  filter(cluster=="PML bodies") %>%
  select(mapobject_id,`11_PML`,`20_SP100`,Nuclei_Intensity_mean_00_EU,cell_cycle) %>%
  group_by(cell_cycle) %>%
  group_by(cell_cycle) %>%
  summarise(cor_PML_EU_all = list(bootstrap_cor_all(`11_PML`,Nuclei_Intensity_mean_00_EU)),
            cor_SP100_EU_all = list(bootstrap_cor_all(`20_SP100`,Nuclei_Intensity_mean_00_EU))) %>%
  unnest(col = c(cor_PML_EU_all,cor_SP100_EU_all),names_sep="_")

correlations_with_EU_all %>%
  pivot_longer(-cell_cycle,values_to="cor") %>%
  mutate(channel = if_else(grepl("PML",name),"PML","SP100")) %>%
  select(-name) %>%
  mutate(cell_cycle = factor(cell_cycle,levels=c("G1","S","G2"))) %>%
  ggplot(aes(x=channel,y=cor,fill=cell_cycle)) +
  geom_violin(position = position_dodge(0.7),width=0.7,col="black",size=0.2) +
  geom_point(data=correlations_with_EU_summary,aes(y=cor),position = position_dodge(0.7),col="black",size=0.2) +
  geom_errorbar(data=correlations_with_EU_summary,aes(ymin=lower,ymax=upper),position = position_dodge(0.7),size=0.3,width=0.5) +
  geom_hline(yintercept=0,size=0.2) +
  scale_y_continuous("Correlation with EU",limits=c(-.5,.5),breaks=c(-.4,0,.4)) +
  scale_fill_brewer(name="Cell cycle\nstage", palette = "Greys") +
  theme_bw(base_size = 7) +
  #ggtitle("PML bodies") +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_text(size=6),
        #axis.text.x = element_text(angle=90,hjust=1,vjust=0.5),
        legend.key.size = unit(2,"mm"),
        legend.position = "right",
        legend.margin = margin(-4,0,-4,-4),
        legend.direction = "vertical")
ggsave(file.path(plot_dir,"PML_SP100_EU_violin_plot.pdf"),width=4.2,height=3,units="cm")


# Example cells ----

# This is very clear in S-phase cells.
# Find an S-phase cell from each bin to make the comparison

manual_selection <- tibble(mapobject_id = c(334228,208743))
manual_selection_ids <- unique(manual_selection$mapobject_id)

pixels <- campa_res %>%
  inner_join(manual_selection) %>%
  filter(!is.na(campa_res_dir)) %>%
  distinct(campa_res_dir,data_dir) %>%
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
  filter(mapobject_id %in% manual_selection_ids) %>%
  mutate(mapobject_id = as.integer(mapobject_id)) %>%
  select(-campa_res_dir,-data_dir) %>%
  data.table::as.data.table()

# get pixel clustering only
campa_pixels <- pixels %>%
  select(-input_data,-channels) %>%
  unnest(leiden_cluster_id) %>%
  filter(mapobject_id %in% manual_selection_ids) %>%
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
  tibble() 

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

# plot just these cells and save plot
segmentation_images <- manual_selection %>%
  left_join(selected_pixels) %>%
  mutate(x_new = if_else(mapobject_id==334228,x,y),
         y_new = if_else(mapobject_id==334228,y,x),
         x=x_new,
         y=y_new) %>%
  select(-x_new,-y_new) %>%
  ggplot(aes(x=x,y=y,fill=cluster_annotation)) +
  geom_tile() +
  coord_fixed() +
  facet_wrap(~mapobject_id) +
  scale_fill_manual(values=getcolor[levels(campa_pixels$cluster_annotation)]) +
  theme_void(base_size = 7) +
  theme(legend.title = element_blank(),
        legend.key.size = unit(1.5,"mm"))
segmentation_images

ggsave(plot = segmentation_images,
       filename = file.path(plot_dir,"EU_selected_segmentation.pdf"),
       width=12,height=7,units="cm")
ggsave_cairo(plot = segmentation_images,
             filename = file.path(plot_dir,"EU_selected_segmentation.png"),
             width=12,height=7,units="cm")

# plot just these cells and save plot
segmentation_images_PML <- manual_selection %>%
  left_join(selected_pixels) %>%
  mutate(x_new = if_else(mapobject_id==334228,x,y),
         y_new = if_else(mapobject_id==334228,y,x),
         x=x_new,
         y=y_new) %>%
  select(-x_new,-y_new) %>%
  mutate(cluster_PML = factor(case_when(cluster_annotation=="PML bodies" ~ "PML",
                                        cluster_annotation %in% nuclear_CSLs ~ "Nucleus",
                                                   TRUE ~ "Cytoplasm"),levels=c("PML","Nucleus","Cytoplasm"))) %>%
  ggplot(aes(x=x,y=y,fill=cluster_PML)) +
  geom_tile() +
  coord_fixed() +
  facet_wrap(~mapobject_id) +
  scale_fill_manual(values=c(getcolor["PML bodies"][[1]],"grey80","grey90")) +
  theme_void(base_size = 7) +
  theme(legend.title = element_blank(),
        legend.key.size = unit(1.5,"mm"))
segmentation_images_PML

ggsave(plot = segmentation_images_PML,
       filename = file.path(plot_dir,"EU_selected_segmentation_PML.pdf"),
       width=12,height=7,units="cm")
ggsave_cairo(plot = segmentation_images_PML,
             filename = file.path(plot_dir,"EU_selected_segmentation_PML.png"),
             width=12,height=7,units="cm")


# plot intensity images
intensity_images <- manual_selection %>%
  select(mapobject_id) %>%
  left_join(selected_pixels) %>%
  ungroup() %>%
  mutate(x_new = if_else(mapobject_id==334228,x,y),
         y_new = if_else(mapobject_id==334228,y,x),
         x=x_new,
         y=y_new) %>%
  select(-x_new,-y_new) %>%
  select(mapobject_id,x,y,`11_PML`) %>%
  mutate_at(vars(matches("\\d{2}_")),~./ quantile(.,0.9995)) %>%
  pivot_longer(cols=matches("\\d{2}_")) %>%
  ggplot(aes(x=x,y=y,fill=value)) +
  geom_tile() +
  coord_fixed() +
  facet_wrap(mapobject_id~name) +
  scale_fill_gradient(low="white",high="black",limits=c(0,1.2),breaks=c(0,0.5,1),oob=scales::squish) +
  theme_void(base_size = 7) +
  theme(legend.title = element_blank(),
        legend.key.size = unit(1.5,"mm"))
intensity_images

ggsave(plot = intensity_images,
       filename = file.path(plot_dir,"EU_selected_PML.pdf"),
       width=7,height=3,units="cm")
ggsave_cairo(plot = intensity_images,
             filename = file.path(plot_dir,"EU_selected_PML.png"),
             width=7,height=3,units="cm")

intensity_images <- manual_selection %>%
  select(mapobject_id) %>%
  left_join(selected_pixels) %>%
  ungroup() %>%
  mutate(x_new = if_else(mapobject_id==334228,x,y),
         y_new = if_else(mapobject_id==334228,y,x),
         x=x_new,
         y=y_new) %>%
  select(-x_new,-y_new) %>%
  select(mapobject_id,x,y,`20_SP100`) %>%
  mutate_at(vars(matches("\\d{2}_")),~./ quantile(.,0.999)) %>%
  pivot_longer(cols=matches("\\d{2}_")) %>%
  ggplot(aes(x=x,y=y,fill=value)) +
  geom_tile() +
  coord_fixed() +
  facet_wrap(mapobject_id~name) +
  scale_fill_gradient(low="white",high="black",limits=c(0,1.2),breaks=c(0,0.5,1),oob=scales::squish) +
  theme_void(base_size = 7) +
  theme(legend.title = element_blank(),
        legend.key.size = unit(1.5,"mm"))
intensity_images

ggsave(plot = intensity_images,
       filename = file.path(plot_dir,"EU_selected_SP100.pdf"),
       width=7,height=3,units="cm")
ggsave_cairo(plot = intensity_images,
             filename = file.path(plot_dir,"EU_selected_SP100.png"),
             width=7,height=3,units="cm")

# Add outlines
dummy_pixels <- manual_selection %>%
  left_join(selected_pixels) %>%
  group_by(mapobject_id) %>%
  select(x,y) %>%
  pivot_longer(cols=c(x,y)) %>%
  group_by(mapobject_id,name) %>%
  summarise(min = min(value),
            max = max(value)) %>%
  pivot_wider(values_from = c(min,max))

dummy_pixels <- bind_rows(select(dummy_pixels,mapobject_id,x=min_x,y=min_y),
                          select(dummy_pixels,mapobject_id,x=max_x,y=max_y))

condvae_pml_outlines <- manual_selection %>%
  left_join(selected_pixels) %>%
  filter(cluster_annotation %in% c("PML bodies")) %>%
  select(mapobject_id,x,y) %>%
  bind_rows(dummy_pixels) %>%
  group_by(mapobject_id) %>%
  nest() %>%
  mutate(sparse = map(data,to_sparse)) %>%
  mutate(sparse = map(sparse,~dilate_sparse(.,width=9))) %>%
  mutate(outlines = map(sparse,~to_outlines_sparse(.,size=1))) %>%
  select(-data,-sparse) %>%
  unnest(cols=outlines) %>%
  mutate(outline_pml=T,y_new = x - 1, x_new = y - 1) %>%
  select(outline_pml,y=y_new,x=x_new)

condvae_all_outlines <- manual_selection %>%
  left_join(selected_pixels) %>%
  select(mapobject_id,x,y) %>%
  bind_rows(dummy_pixels) %>%
  group_by(mapobject_id) %>%
  nest() %>%
  mutate(sparse = map(data,to_sparse)) %>%
  mutate(sparse = map(sparse,~dilate_sparse(.,width=3))) %>%
  mutate(outlines = map(sparse,~to_outlines_sparse(.,size=3))) %>%
  select(-data,-sparse) %>%
  unnest(cols=outlines) %>%
  mutate(outline_all=T,y_new = x - 1, x_new = y - 1) %>%
  select(outline_all,y=y_new,x=x_new)

# plot intensity images
intensity_images <- manual_selection %>%
  left_join(selected_pixels) %>%
  left_join(condvae_pml_outlines) %>%
  left_join(condvae_all_outlines) %>%
  ungroup() %>%
  select(-`00_EU`) %>%
  mutate(x_new = if_else(mapobject_id==334228,x,y),
         y_new = if_else(mapobject_id==334228,y,x),
         x=x_new,
         y=y_new) %>%
  select(-x_new,-y_new) %>%
  mutate_at(vars(matches("\\d{2}_")),~./ quantile(.,0.9995)) %>% {
    ggplot(data=.,aes(x=x,y=y,fill=`11_PML`)) +
      geom_tile() +
      geom_tile(data=filter(.,outline_pml==T),fill=getcolor["PML bodies"][[1]],alpha=0.5) +
      geom_tile(data=filter(.,outline_all==T),fill="grey30") } +
  coord_equal(ratio = 1) +
  facet_wrap(~mapobject_id) +
  scale_fill_gradient(low="white",high="black",limits=c(0,1.2),breaks=c(0,0.5,1),oob=scales::squish) +
  theme_void(base_size = 7) +
  theme(legend.title = element_blank(),
        legend.key.size = unit(1.5,"mm"))

intensity_images
ggsave(plot = intensity_images,
       filename = file.path(plot_dir,"EU_selected_PML_outlines.pdf"),
       width=7,height=3,units="cm")
ggsave_cairo(plot = intensity_images,
             filename = file.path(plot_dir,"EU_selected_PML_outlines.png"),
             width=7,height=3,units="cm",dpi = 600)


# plot intensity images
intensity_images <- manual_selection %>%
  left_join(selected_pixels) %>%
  left_join(condvae_pml_outlines) %>%
  left_join(condvae_all_outlines) %>%
  ungroup() %>%
  select(-`00_EU`) %>%
  mutate(x_new = if_else(mapobject_id==334228,x,y),
         y_new = if_else(mapobject_id==334228,y,x),
         x=x_new,
         y=y_new) %>%
  select(-x_new,-y_new) %>%
  mutate_at(vars(matches("\\d{2}_")),~./ quantile(.,0.9995)) %>% {
    ggplot(data=.,aes(x=x,y=y,fill=`20_SP100`)) +
      geom_tile() +
      geom_tile(data=filter(.,outline_pml==T),fill=getcolor["PML bodies"][[1]],alpha=0.5) +
      geom_tile(data=filter(.,outline_all==T),fill="grey30") } +
  coord_equal(ratio = 1) +
  facet_wrap(~mapobject_id) +
  scale_fill_gradient(low="white",high="black",limits=c(0,1.2),breaks=c(0,0.5,1),oob=scales::squish) +
  theme_void(base_size = 7) +
  theme(legend.title = element_blank(),
        legend.key.size = unit(1.5,"mm"))

intensity_images
ggsave(plot = intensity_images,
       filename = file.path(plot_dir,"EU_selected_SP100_outlines.pdf"),
       width=7,height=3,units="cm")
ggsave_cairo(plot = intensity_images,
             filename = file.path(plot_dir,"EU_selected_SP100_outlines.png"),
             width=7,height=3,units="cm",dpi = 600)


