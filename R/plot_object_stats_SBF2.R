# clear the workspace
rm(list=ls())

# load required packages
library(tidyverse)
library(patchwork)
library(ComplexHeatmap)

# setup script-specific parameters
cell_type <- "HeLa"
experiment_name <- "VAE_SBF2/CondVAE_siRNA-CC"
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
plot_dir <- file.path(campa_ana$constants$SOURCE_DIR,"figures","object_stats")
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
}

# read cluster annotation
annotation <- read_csv(file.path(campa_ana$constants$SOURCE_DIR,"R","annotation_VAE_SBF2.csv"))

# make a color look-up table (annotations)
getcolor <- distinct(annotation,cluster,color) %>% pull(color)
names(getcolor) <- distinct(annotation,cluster,color) %>% pull(cluster)

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
  inner_join(data_dirs,by="well_name")

# load campa results
campa_res <- selected_wells %>%
  mutate(model_res = map(campa_res_dir,~read_csv(file.path(.,"features_seed3_annotation.csv"),show_col_types = FALSE) %>% select(-1))) %>%
  select(-well_name,-siRNA) %>%
  unnest(model_res) %>%
  mutate(siRNA = factor(siRNA,levels=c("scrambled","SBF2")))

# load spatial stats
object_stats <- selected_wells %>%
  mutate(object_stats = map(campa_res_dir,~read_csv(file.path(.,"export","object_stats_features_seed3_annotation.csv"),show_col_types = FALSE) %>% select(-1))) %>%
  select(-well_name) %>%
  unnest(object_stats) %>%
  pivot_longer(cols=matches("\\|"),names_to = c("statistic","cluster"),names_sep="\\|") %>%
  mutate(siRNA = factor(siRNA,levels=c("scrambled","SBF2")))

nuclear_CSLs <- c(
  "Cajal bodies",             
  "Nuclear periphery",        
  "Nuclear speckles",         
  "Nucleolus",                
  "Nucleoplasm",
  "PML bodies")

campa_res_nucleus_or_cytoplasm <- campa_res %>%
  filter(cluster!="all") %>%
  mutate(nucleus_or_cytoplasm = if_else(cluster %in% nuclear_CSLs,"Nucleus (combined)","Cytoplasm (combined)")) %>%
  mutate_at(vars(matches("\\d{2}_")),~.*size) %>%
  group_by(mapobject_id,well_name,siRNA,TR,cell_cycle,nucleus_or_cytoplasm) %>%
  summarise_at(vars(matches("\\d{2}_|size")),~sum(.,na.rm=T)) %>%
  mutate_at(vars(matches("\\d{2}_")),~./size) %>%
  ungroup() %>%
  rename(cluster = nucleus_or_cytoplasm)

selected_sizes_to_plot <- campa_res %>%
  mutate(siRNA = factor(siRNA,levels=c("scrambled","SBF2"))) %>%
  select(well_name,siRNA,mapobject_id,cluster,size) %>%
  left_join(filter(campa_res_nucleus_or_cytoplasm,cluster=="Nucleus (combined)") %>%
              select(nucleus_size=size,mapobject_id)) %>%
  filter(cluster %in% c("Cajal bodies","Nucleolus")) %>%
  mutate(fraction_of_nucleus = size/nucleus_size)

# Compute the confidence interval of the mean co-occurrence score using bootsrapping 
bootstrap_mean_ci <- function(x) {
  out <- boot::boot(data = x,
                    R = 100,
                    statistic = function(dat, idx) mean(dat[idx], na.rm = TRUE))
  ci <- boot::boot.ci(out, type = "perc",conf = 0.95)
  return(tibble(lower = ci$percent[4],upper = ci$percent[5]))
}

boxplot_theme <- theme_bw(base_size = 7) +
  theme(legend.position = "none",
        panel.grid=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=6))

# Nucleolus
nucleolus_size_absolute <- campa_res %>%
  filter(cluster == "Nucleolus") %>%
  ggplot(aes(y=size*pixel_size*pixel_size,fill=siRNA)) + 
  geom_boxplot(outlier.shape = NA,size=0.3) +
  scale_fill_manual(values=c("#999999",cbPalette[3])) +
  boxplot_theme +
  coord_cartesian(ylim=c(0,100)) +
  scale_y_continuous(name='Total area',breaks=scales::pretty_breaks(n=3))
nucleolus_size_absolute

nucleolus_size_fraction <- selected_sizes_to_plot %>%
  filter(cluster=="Nucleolus") %>%
  ggplot(aes(y=fraction_of_nucleus,fill=siRNA)) + 
  geom_boxplot(outlier.shape = NA,size=0.3) +
  scale_fill_manual(values=c("#999999",cbPalette[3])) +
  boxplot_theme +
  coord_cartesian(ylim=c(0,0.4)) +
  scale_y_continuous(name="Area (% of nucleus)",breaks=scales::pretty_breaks(n=3),labels = ~scales::percent(.,accuracy = 1))
nucleolus_size_fraction

NCL_abundance_data <- campa_res_nucleus_or_cytoplasm %>%
  mutate(siRNA = factor(siRNA,levels=c("scrambled","SBF2"))) %>%
  filter(cluster=="Nucleus (combined)") %>%
  mutate(sum_NCL = size*`21_NCL`) 

NCL_abundance <- NCL_abundance_data %>%
  filter(siRNA=="scrambled") %>%
  summarise(normaliser = median(sum_NCL)) %>%
  bind_cols(NCL_abundance_data) %>%
  ggplot(aes(y=sum_NCL/normaliser,fill=siRNA)) + 
  geom_boxplot(outlier.shape = NA,size=0.3) +
  scale_fill_manual(values=c("#999999",cbPalette[3])) +
  boxplot_theme +
  coord_cartesian(ylim=c(0,5)) +
  scale_y_continuous(name="NCL (Nucleus)",breaks=scales::pretty_breaks(n=3))
NCL_abundance

NCL_abundance_in_nucleolus_data <- campa_res %>%
  mutate(siRNA = factor(siRNA,levels=c("scrambled","SBF2"))) %>%
  filter(cluster=="Nucleolus") %>%
  mutate(sum_NCL = size*`21_NCL`) 

NCL_abundance_in_nucleolus <- NCL_abundance_in_nucleolus_data %>%
  filter(siRNA=="scrambled") %>%
  summarise(normaliser = median(sum_NCL)) %>%
  bind_cols(NCL_abundance_in_nucleolus_data) %>%
  ggplot(aes(y=sum_NCL/normaliser,fill=siRNA)) + 
  geom_boxplot(outlier.shape = NA,size=0.3) +
  scale_fill_manual(values=c("#999999",cbPalette[3])) +
  boxplot_theme +
  coord_cartesian(ylim=c(0,5)) +
  scale_y_continuous(name="NCL (Nucleolus)",breaks=scales::pretty_breaks(n=3))
NCL_abundance_in_nucleolus

# individual 
nucleolus_count_thresholded <- object_stats %>%
  filter(cluster == "Nucleolus" & statistic=="area_thresholded_count") %>% 
  ggplot(aes(y=value,fill=siRNA)) + 
  geom_boxplot(outlier.shape = NA,size=0.3) +
  scale_fill_manual(values=c("#999999",cbPalette[3])) +
  boxplot_theme +
  coord_cartesian(ylim=c(0,5)) +
  scale_y_continuous(name="Nucleoli per cell",breaks=scales::pretty_breaks(n=3))
nucleolus_count_thresholded

nucleolus_area_thresholded <- object_stats %>%
  filter(cluster == "Nucleolus" & statistic=="area_thresholded_median") %>% 
  ggplot(aes(y=value*pixel_size*pixel_size,fill=siRNA)) + 
  geom_boxplot(outlier.shape = NA,size=0.3) +
  scale_fill_manual(values=c("#999999",cbPalette[3])) +
  boxplot_theme +
  coord_cartesian(ylim=c(0,75)) +
  scale_y_continuous(name='Individual\nnucleoli area',breaks=scales::pretty_breaks(n=3))
nucleolus_area_thresholded

nucleolus_count_thresholded + nucleolus_area_thresholded + nucleolus_size_fraction + nucleolus_size_absolute + plot_layout(ncol=2,nrow=2) + theme(plot.margin = margin(0,0,0,0))
ggsave(file.path(plot_dir,"NCL_object_stats_main_combined.pdf"),
       width=3.8,height=4.5,units="cm")

NCL_abundance + NCL_abundance_in_nucleolus + plot_layout(ncol=2,nrow=1) + theme(plot.margin = margin(0,0,0,0))
ggsave(file.path(plot_dir,"NCL_object_stats_si_combined.pdf"),
       width=3.8,height=2.5,units="cm")


# Cajal bodies
cajal_size_absolute <- campa_res %>%
  filter(cluster == "Cajal bodies") %>%
  ggplot(aes(y=size*pixel_size*pixel_size,fill=siRNA)) + 
  geom_boxplot(outlier.shape = NA,size=0.3) +
  scale_fill_manual(values=c("#999999",cbPalette[3])) +
  boxplot_theme +
  coord_cartesian(ylim=c(0,50)) +
  scale_y_continuous(name='Total area',breaks=scales::pretty_breaks(n=3))
cajal_size_absolute

cajal_size_fraction <- selected_sizes_to_plot %>%
  filter(cluster=="Cajal bodies") %>%
  ggplot(aes(y=fraction_of_nucleus,fill=siRNA)) + 
  geom_boxplot(outlier.shape = NA,size=0.3) +
  scale_fill_manual(values=c("#999999",cbPalette[3])) +
  boxplot_theme +
  coord_cartesian(ylim=c(0,0.1)) +
  scale_y_continuous(name="Area (% of nucleus)", breaks=scales::pretty_breaks(n=3),labels = ~scales::percent(.,accuracy = 1))
cajal_size_fraction

COIL_abundance_data <- campa_res_nucleus_or_cytoplasm %>%
  mutate(siRNA = factor(siRNA,levels=c("scrambled","SBF2"))) %>%
  filter(cluster=="Nucleus (combined)") %>%
  mutate(sum_COIL = size*`21_COIL`) 

COIL_abundance <- COIL_abundance_data %>%
  filter(siRNA=="scrambled") %>%
  summarise(normaliser = median(sum_COIL)) %>%
  bind_cols(COIL_abundance_data) %>%
  ggplot(aes(y=sum_COIL/normaliser,fill=siRNA)) + 
  geom_boxplot(outlier.shape = NA,size=0.3) +
  scale_fill_manual(values=c("#999999",cbPalette[3])) +
  boxplot_theme +
  coord_cartesian(ylim=c(0,10)) +
  scale_y_continuous(name="COIL (Nucleus)",breaks=scales::pretty_breaks(n=3))
COIL_abundance

COIL_abundance_in_cajal_data <- campa_res %>%
  mutate(siRNA = factor(siRNA,levels=c("scrambled","SBF2"))) %>%
  filter(cluster=="Cajal bodies") %>%
  mutate(sum_COIL = size*`21_COIL`) 

COIL_abundance_in_cajal <- COIL_abundance_in_cajal_data %>%
  filter(siRNA=="scrambled") %>%
  summarise(normaliser = median(sum_COIL)) %>%
  bind_cols(COIL_abundance_in_cajal_data) %>%
  ggplot(aes(y=sum_COIL/normaliser,fill=siRNA)) + 
  geom_boxplot(outlier.shape = NA,size=0.3) +
  scale_fill_manual(values=c("#999999",cbPalette[3])) +
  boxplot_theme +
  coord_cartesian(ylim=c(0,10)) +
  scale_y_continuous(name="COIL (Cajal bodies)",breaks=scales::pretty_breaks(n=3))
COIL_abundance_in_cajal

# individual 
cajal_count_thresholded <- object_stats %>%
  filter(cluster == "Cajal bodies" & statistic=="area_thresholded_count") %>% 
  ggplot(aes(y=value,fill=siRNA)) + 
  geom_boxplot(outlier.shape = NA,size=0.3) +
  scale_fill_manual(values=c("#999999",cbPalette[3])) +
  boxplot_theme +
  coord_cartesian(ylim=c(0,10)) +
  scale_y_continuous(name="Cajal bodies per cell",breaks=scales::pretty_breaks(n=3))
cajal_count_thresholded

cajal_area_thresholded <- object_stats %>%
  filter(cluster == "Cajal bodies" & statistic=="area_thresholded_median") %>% 
  ggplot(aes(y=value*pixel_size*pixel_size,fill=siRNA)) + 
  geom_boxplot(outlier.shape = NA,size=0.3) +
  scale_fill_manual(values=c("#999999",cbPalette[3])) +
  boxplot_theme +
  coord_cartesian(ylim=c(0,20)) +
  scale_y_continuous(name='Individual\nCajal-body area',breaks=scales::pretty_breaks(n=3))
cajal_area_thresholded

cajal_count_thresholded + cajal_area_thresholded + cajal_size_fraction + cajal_size_absolute + plot_layout(ncol=2,nrow=2) + theme(plot.margin = margin(0,0,0,0))
ggsave(file.path(plot_dir,"COIL_object_stats_main_combined.pdf"),
       width=3.8,height=4.5,units="cm")

COIL_abundance + COIL_abundance_in_cajal + plot_layout(ncol=2,nrow=1) + theme(plot.margin = margin(0,0,0,0))
ggsave(file.path(plot_dir,"COIL_object_stats_si_combined.pdf"),
       width=3.8,height=2.5,units="cm")

# P-bodies
selected_sizes_to_plot_cytoplasm <- campa_res %>%
  mutate(siRNA = factor(siRNA,levels=c("scrambled","SBF2"))) %>%
  select(well_name,siRNA,mapobject_id,cluster,size) %>%
  left_join(filter(campa_res_nucleus_or_cytoplasm,cluster=="Cytoplasm (combined)") %>%
              select(cytoplasm_size=size,mapobject_id)) %>%
  filter(cluster %in% c("ER (perinuclear)","P-bodies")) %>%
  mutate(fraction_of_cytoplasm = size/cytoplasm_size)

p_body_size_absolute <- campa_res %>%
  filter(cluster == "P-bodies") %>%
  ggplot(aes(y=size*pixel_size*pixel_size,fill=siRNA)) + 
  geom_boxplot(outlier.shape = NA,size=0.3) +
  scale_fill_manual(values=c("#999999",cbPalette[3])) +
  boxplot_theme +
  coord_cartesian(ylim=c(0,120)) +
  scale_y_continuous(name='Total area',breaks=scales::pretty_breaks(n=3))
p_body_size_absolute

p_body_size_fraction <- selected_sizes_to_plot_cytoplasm %>%
  filter(cluster=="P-bodies") %>%
  ggplot(aes(y=fraction_of_cytoplasm,fill=siRNA)) + 
  geom_boxplot(outlier.shape = NA,size=0.3) +
  scale_fill_manual(values=c("#999999",cbPalette[3])) +
  boxplot_theme +
  coord_cartesian(ylim=c(0,0.15)) +
  scale_y_continuous(name="Area (% of cytoplasm)",breaks=scales::pretty_breaks(n=3),labels = ~scales::percent(.,accuracy = 1))
p_body_size_fraction

DDX6_abundance_data <- campa_res_nucleus_or_cytoplasm %>%
  mutate(siRNA = factor(siRNA,levels=c("scrambled","SBF2"))) %>%
  filter(cluster=="Nucleus (combined)") %>%
  mutate(sum_DDX6 = size*`22_DDX6`) 

DDX6_abundance <- DDX6_abundance_data %>%
  filter(siRNA=="scrambled") %>%
  summarise(normaliser = median(sum_DDX6)) %>%
  bind_cols(DDX6_abundance_data) %>%
  ggplot(aes(y=sum_DDX6/normaliser,fill=siRNA)) + 
  geom_boxplot(outlier.shape = NA,size=0.3) +
  scale_fill_manual(values=c("#999999",cbPalette[3])) +
  boxplot_theme +
  coord_cartesian(ylim=c(0,4)) +
  scale_y_continuous(name="DDX6 (Cytoplasm)", breaks=scales::pretty_breaks(n=3))
DDX6_abundance

DDX6_abundance_in_p_body_data <- campa_res %>%
  mutate(siRNA = factor(siRNA,levels=c("scrambled","SBF2"))) %>%
  filter(cluster=="P-bodies") %>%
  mutate(sum_DDX6 = size*`22_DDX6`) 

DDX6_abundance_in_p_body <- DDX6_abundance_in_p_body_data %>%
  filter(siRNA=="scrambled") %>%
  summarise(normaliser = median(sum_DDX6)) %>%
  bind_cols(DDX6_abundance_in_p_body_data) %>%
  ggplot(aes(y=sum_DDX6/normaliser,fill=siRNA)) + 
  geom_boxplot(outlier.shape = NA,size=0.3) +
  scale_fill_manual(values=c("#999999",cbPalette[3])) +
  boxplot_theme +
  coord_cartesian(ylim=c(0,4)) +
  scale_y_continuous(name="DDX6 (P-bodies)",breaks=scales::pretty_breaks(n=3))
DDX6_abundance_in_p_body

# individual 
p_body_count_thresholded <- object_stats %>%
  filter(cluster == "P-bodies" & statistic=="area_thresholded_count") %>% 
  ggplot(aes(y=value,fill=siRNA)) + 
  geom_boxplot(outlier.shape = NA,size=0.3) +
  scale_fill_manual(values=c("#999999",cbPalette[3])) +
  boxplot_theme +
  coord_cartesian(ylim=c(0,120)) +
  scale_y_continuous(name="P-bodies per cell",breaks=scales::pretty_breaks(n=3))
p_body_count_thresholded

p_body_area_thresholded <- object_stats %>%
  filter(cluster == "P-bodies" & statistic=="area_thresholded_median") %>% 
  ggplot(aes(y=value*pixel_size*pixel_size,fill=siRNA)) + 
  geom_boxplot(outlier.shape = NA,size=0.3) +
  scale_fill_manual(values=c("#999999",cbPalette[3])) +
  boxplot_theme +
  coord_cartesian(ylim=c(0,2.1)) +
  scale_y_continuous(name='Individual\nP-body area',breaks=scales::pretty_breaks(n=3))
p_body_area_thresholded
 
p_body_count_thresholded + p_body_area_thresholded + p_body_size_fraction + p_body_size_absolute + plot_layout(ncol=2,nrow=2) + theme(plot.margin = margin(0,0,0,0))
ggsave(file.path(plot_dir,"DDX6_object_stats_main_combined.pdf"),
       width=3.8,height=4.5,units="cm")

DDX6_abundance + DDX6_abundance_in_p_body + plot_layout(ncol=2,nrow=1) + theme(plot.margin = margin(0,0,0,0))
ggsave(file.path(plot_dir,"DDX6_object_stats_si_combined.pdf"),
       width=3.8,height=2.5,units="cm")

# plot P-bodies by cell size bins

cell_features <- read_csv(file=file.path(campa_ana$constants$SOURCE_DIR,"additional_data","HeLa_cell_features.csv"))

cell_features %>%
  left_join(select(campa_res,mapobject_id,siRNA)) %>%
  ggplot(aes(x=Intensity_sum_22_SE,group=well_name,col=siRNA)) + geom_density()

# bin by protein content
cell_features <- cell_features %>%
  ungroup() %>%
  mutate(size_bin = cut_interval(Intensity_sum_22_SE,n=30,labels=F))

p_body_count_data <- object_stats %>%
  filter(cluster == "P-bodies") %>% 
  filter(statistic=="area_thresholded_count") %>%
  left_join(cell_features) %>%
  mutate(siRNA = factor(siRNA,levels=c("scrambled","SBF2")))

p_counts <- p_body_count_data %>%
  group_by(siRNA,size_bin) %>%
  count() %>%
  filter(n>1) %>%
  inner_join(p_body_count_data) %>%
  group_by(siRNA,size_bin) %>%
  summarise(n = n(),
            mean = mean(value),
            mean_ci = list(bootstrap_mean_ci(value))) %>%
  filter(n > 5) %>%
  left_join( distinct(p_body_count_data,mapobject_id,siRNA) %>% count(siRNA,name="total_n")) %>%
  unnest_wider(mean_ci) %>%
  ungroup() %>%
  pivot_wider(names_from = siRNA,values_from = n) %>%
  mutate(across(SBF2:scrambled,~replace_na(.x,0))) %>%
  pivot_longer(SBF2:scrambled,names_to="siRNA",values_to = "n") %>%
  mutate(siRNA = factor(siRNA,levels=c("scrambled","SBF2"))) %>%
  ggplot(aes(x=size_bin-1,y=n/total_n,fill=siRNA)) +
  geom_col(position = "dodge") +
  scale_y_continuous(breaks = c(0,0.1)) +
  scale_fill_manual(values=c("#999999",cbPalette[3])) +
  coord_cartesian(ylim=c(0,NA),xlim=c(-1,20),expand = c(0,0)) +
  ylab("Fraction\nof cells") +
  theme_bw(base_size = 7) +
  theme(panel.border = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.y = element_line(size=0.2),
        axis.line.x = element_line(size=0.2),
        panel.spacing = unit(1,"mm"),
        panel.grid = element_blank(),
        legend.position = "none",
        axis.title=element_text(size=6))
p_counts

p_trends_count <- p_body_count_data %>%
  group_by(siRNA,size_bin) %>%
  count() %>%
  filter(n>1) %>%
  inner_join(p_body_count_data) %>%
  group_by(siRNA,size_bin) %>%
  summarise(n = n(),
            mean = mean(value),
            mean_ci = list(bootstrap_mean_ci(value))) %>%
  filter(n > 5) %>%
  unnest_wider(mean_ci) %>%
  ungroup() %>%
  mutate(siRNA = factor(siRNA,levels=c("scrambled","SBF2"))) %>%
  ggplot(aes(x=size_bin-1,y=mean,col=siRNA)) +
  geom_errorbar(aes(ymin=lower,ymax=upper),width=0.2,size=0.2) +
  geom_smooth(data=p_body_count_data %>%
                group_by(siRNA) %>%
                filter(Intensity_sum_22_SE < quantile(Intensity_sum_22_SE,0.995) &
                         Intensity_sum_22_SE > quantile(Intensity_sum_22_SE,0.005)),
              aes(y=value,fill=siRNA),size=0.2,alpha=0.2) +
  geom_point(pch=21,fill="white",size=0.5) +
  expand_limits(x=0,y=0) +
  scale_color_manual(values=c("#999999",cbPalette[3])) +
  scale_fill_manual(values=c("#999999",cbPalette[3])) +
  scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
  coord_cartesian(ylim=c(0,NA),xlim=c(-1,20),expand = c(0,0)) +
  ylab("P-bodies\nper cell") +
  xlab("Cell size (total protein content) bin") +
  theme_bw(base_size = 7) +
  theme(panel.spacing = unit(1,"mm"),
        panel.grid = element_blank(),
        legend.title=element_text(size=6),
        legend.key.size = unit(2,"mm"),
        axis.title=element_text(size=6))

combined_plot <- p_counts / p_trends_count + plot_layout(heights = c(1,4))
combined_plot
ggsave(plot = combined_plot, filename = file.path(plot_dir,"P-bodies_count_trend_SE.pdf"),width=6,height=3.8,units="cm")
ggsave_cairo(plot = combined_plot, filename = file.path(plot_dir,"P-bodies_count_trend_SE.png"),width=6,height=3.8,units="cm",dpi=600)

p_body_size_data <- object_stats %>%
  filter(cluster == "P-bodies") %>% 
  filter(statistic=="area_thresholded_median") %>%
  mutate(value=value*pixel_size*pixel_size) %>%
  left_join(cell_features) %>%
  mutate(siRNA = factor(siRNA,levels=c("scrambled","SBF2")))

p_trends_size <- p_body_size_data %>%
  group_by(siRNA,size_bin) %>%
  count() %>%
  filter(n>1) %>%
  inner_join(p_body_size_data) %>%
  group_by(siRNA,size_bin) %>%
  summarise(n = n(),
            mean = mean(value),
            mean_ci = list(bootstrap_mean_ci(value))) %>%
  filter(n > 5) %>%
  unnest_wider(mean_ci) %>%
  ungroup() %>%
  mutate(siRNA = factor(siRNA,levels=c("scrambled","SBF2"))) %>%
  ggplot(aes(x=size_bin-1,y=mean,col=siRNA)) +
  geom_errorbar(aes(ymin=lower,ymax=upper),width=0.2,size=0.2) +
  geom_smooth(data=p_body_size_data %>%
                group_by(siRNA) %>%
                filter(Intensity_sum_22_SE < quantile(Intensity_sum_22_SE,0.995) &
                         Intensity_sum_22_SE > quantile(Intensity_sum_22_SE,0.005)),
              aes(y=value,fill=siRNA),size=0.2,alpha=0.2) +
  geom_point(pch=21,fill="white",size=0.5) +
  expand_limits(x=0,y=0) +
  scale_color_manual(values=c("#999999",cbPalette[3])) +
  scale_fill_manual(values=c("#999999",cbPalette[3])) +
  scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
  coord_cartesian(ylim=c(0,NA),xlim=c(-1,20),expand = c(0,0)) +
  #ylab(bquote('P-body area('*µm^2*')')) +
  ylab("Individual\nP-body area") +
  xlab("Cell size (total protein content) bin") +
  theme_bw(base_size = 7) +
  theme(panel.spacing = unit(1,"mm"),
        panel.grid = element_blank(),
        legend.title=element_text(size=6),
        legend.key.size = unit(2,"mm"),
        axis.title=element_text(size=6))
p_trends_size

combined_plot <- p_counts / p_trends_size + plot_layout(heights = c(1,4))
combined_plot
ggsave(plot = combined_plot, filename = file.path(plot_dir,"P-bodies_ares_trend_SE.pdf"),width=6,height=3.8,units="cm")
ggsave_cairo(plot = combined_plot, filename = file.path(plot_dir,"P-bodies_area_trend_SE.png"),width=6,height=3.8,units="cm",dpi=600)

combined_plot_all <- p_counts / (p_trends_count + theme(axis.text.x=element_blank(),axis.title.x = element_blank())) / p_trends_size +
  plot_layout(heights = c(1,4,4),guides = "collect")
combined_plot_all
ggsave(plot = combined_plot_all, filename = file.path(plot_dir,"P-bodies_count_area_trend_SE.pdf"),width=6,height=5.5,units="cm")
ggsave_cairo(plot = combined_plot_all, filename = file.path(plot_dir,"P-bodies_count_area_trend_SE.png"),width=6,height=5.5,units="cm",dpi=600)
