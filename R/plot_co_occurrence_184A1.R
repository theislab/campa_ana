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
plot_dir <- file.path(campa_ana$constants$SOURCE_DIR,"figures","co_occurrence")
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
}

# read cluster annotation
annotation <- read_csv(file.path(campa_ana$constants$SOURCE_DIR,"R","annotation_VAE_all.csv"))
# read cluster annotation
perturbations <- read_csv(file.path(campa_ana$constants$SOURCE_DIR,"R","perturbation_colors_names.csv"))

# make a color look-up table (annotations)
getcolor <- distinct(annotation,cluster_annotation,color) %>% pull(color)
names(getcolor) <- distinct(annotation,cluster_annotation,color) %>% pull(cluster_annotation)

# make a color look-up table (perturbations)
getcolor_perturbation <- distinct(perturbations,treatment_name,color) %>% pull(color)
names(getcolor_perturbation) <- distinct(perturbations,treatment_name,color) %>% pull(treatment_name)

# read metadata files
channels_metadata <- read_csv(file.path(DATA_DIR,"channels_metadata.csv")) %>% 
  select(-1)
wells_metadata <- read_csv(file.path(DATA_DIR,"wells_metadata.csv")) %>% 
  select(-1) %>%
  mutate(treatment=if_else(perturbation_duration %in% c("normal","DMSO-120","DMSO-720"),
                           "Unperturbed",
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
  inner_join(campa_res_dirs,by="well_name") %>%
  inner_join(data_dirs,by="well_name")

# load co-occurrence scores
co_occurrences <- selected_wells %>%
  mutate(co_occurrence_file = map(campa_res_dir,~list.files(path=file.path(.,"export")))) %>%
  unnest(co_occurrence_file) %>%
  filter(grepl("co_occ",co_occurrence_file)) %>%
  filter(grepl("features_annotation.csv",co_occurrence_file)) %>%
  mutate(co_occ_string = str_remove(string = co_occurrence_file, pattern = "co_occurrence_")) %>%
  mutate(co_occ_string = str_remove(string = co_occ_string, pattern = "_features_annotation.csv")) %>%
  separate(co_occ_string,into=c("from","to"),sep="_") %>%
  mutate(co_occurrence_res = map2(campa_res_dir,co_occurrence_file,~read_csv(file.path(.x,"export",.y),show_col_types = FALSE) %>% select(-1))) %>%
  select(-well_name) %>%
  unnest(co_occurrence_res)

# load spatial stats
object_stats <- selected_wells %>%
  mutate(object_stats = map(campa_res_dir,~read_csv(file.path(.,"export","object_stats_features_annotation.csv"),show_col_types = FALSE) %>% select(-1))) %>%
  select(-well_name) %>%
  unnest(object_stats)

# Compute the confidence interval of the mean co-occurrence score using bootsrapping 
bootstrap_mean_ci <- function(x) {
  out <- boot::boot(data = x,
                    R = 100,
                    statistic = function(dat, idx) mean(dat[idx], na.rm = TRUE))
  ci <- boot::boot.ci(out, type = "perc",conf = 0.95)
  return(tibble(lower = ci$percent[4],upper = ci$percent[5]))
}

intervals <- names(co_occurrences)[grepl(pattern = "^\\d{1}",names(co_occurrences))]

co_occurrences <- co_occurrences %>%
  pivot_longer(cols=one_of(intervals),names_to="distance_interval",values_to="co_occurrence") %>%
  separate(distance_interval,into=c("min_distance","max_distance"),sep="-") %>%
  mutate(across(min_distance:max_distance,as.numeric))

co_occurrences_mean_by_well <- co_occurrences %>%
  # convert pixels to microns
  mutate(min_distance_um = pixel_size*min_distance, max_distance_um = pixel_size*max_distance) %>%
  # remove NA and zero co occurrences
  filter(!is.na(co_occurrence) & co_occurrence > 0) %>%
  # log2 transform
  mutate(log2_co_occurrence = log2(co_occurrence)) %>%
  # get mean and confidence interval for mean using bootstrapping (percentile)
  group_by(from,to,min_distance_um,max_distance_um,treatment,well_name) %>%
  summarise(mean_co_occurrence = mean(co_occurrence),
            mean_log2_co_occurrence = mean(log2_co_occurrence),
            mean_co_occurrence_ci = list(bootstrap_mean_ci(co_occurrence)),
            mean_log2_co_occurrence_ci = list(bootstrap_mean_ci(log2_co_occurrence)) ) %>%
  ungroup() %>%
  unnest(cols=c(mean_co_occurrence_ci,mean_log2_co_occurrence_ci),names_sep="_")

# plot pairwise (unperturbed versus DMSO)
co_occurrences_mean_by_well %>%
  filter(treatment == "Unperturbed") %>%
  select(-treatment) %>%
  left_join(select(wells_metadata,well_name,treatment,perturbation_duration) %>% 
              mutate(treatment=if_else(perturbation_duration=="normal","Untreated",treatment),
                     treatment=if_else(perturbation_duration=="DMSO-120" | perturbation_duration=="DMSO-720","DMSO",treatment))) %>%
  mutate(treatment=factor(treatment)) %>%
  filter(from!="Extra-nuclear" & to !="Extra-nuclear") %>%
  ggplot(aes(x=min_distance_um,y=mean_log2_co_occurrence,col=treatment,fill=treatment,group=well_name)) + 
  geom_ribbon(aes(ymin=mean_log2_co_occurrence_ci_lower,
                  ymax=mean_log2_co_occurrence_ci_upper),alpha=0.05,col=NA) +
  geom_line(size=0.3) +
  scale_y_continuous(name="Log2 spatial co-occurrence") +
  scale_x_log10(name="Distance (µm)",breaks=c(1,10),limits=c(.3,10)) +
  facet_grid(from~to) +
  annotation_logticks(sides = "b",size=0.2,
                      short = unit(0.05, "cm"),
                      mid = unit(0.1, "cm"),
                      long = unit(0.15, "cm")) +
  geom_hline(yintercept = 0,size=0.3) +
  scale_fill_manual(values = c("darkred","#999999")) +
  scale_color_manual(values = c("darkred","#999999")) +
  theme_bw(base_size = 7) +
  theme(legend.title = element_blank(),
        legend.key.size = unit(3,"mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0,"mm"),
        panel.grid = element_blank(),
        panel.spacing = unit(1,"mm"))
ggsave(file.path(plot_dir,"184A1_co_occurrence_Controls_all.pdf"),width=12,height=12.5,units="cm")

# summarise all co-occurrences (averaged by treatment not well)

co_occurrences_mean <- co_occurrences %>%
  # convert pixels to microns
  mutate(min_distance_um = pixel_size*min_distance, max_distance_um = pixel_size*max_distance) %>%
  # remove NA and zero co occurrences
  filter(!is.na(co_occurrence) & co_occurrence > 0) %>%
  # log2 transform
  mutate(log2_co_occurrence = log2(co_occurrence)) %>%
  # get mean and confidence interval for mean using bootstrapping (percentile)
  group_by(from,to,min_distance_um,max_distance_um,treatment) %>%
  summarise(mean_co_occurrence = mean(co_occurrence),
            mean_log2_co_occurrence = mean(log2_co_occurrence),
            mean_co_occurrence_ci = list(bootstrap_mean_ci(co_occurrence)),
            mean_log2_co_occurrence_ci = list(bootstrap_mean_ci(log2_co_occurrence)) ) %>%
  ungroup() %>%
  unnest(cols=c(mean_co_occurrence_ci,mean_log2_co_occurrence_ci),names_sep="_")

# Meayamycin ----

# plot pairwise
to_plot <- co_occurrences_mean %>%
  # focus on Meayamycin and Unperturbed
  filter(treatment %in% c("Meayamycin (12.5h)","Unperturbed")) %>%
  filter(from!="Extra-nuclear" & to !="Extra-nuclear") %>%
  mutate(treatment=factor(treatment))

to_plot %>%
  filter(from!="Extra-nuclear" & to !="Extra-nuclear") %>%
  ggplot(aes(x=min_distance_um,y=mean_log2_co_occurrence,col=treatment,fill=treatment)) + 
  geom_ribbon(aes(ymin=mean_log2_co_occurrence_ci_lower,
                  ymax=mean_log2_co_occurrence_ci_upper),alpha=0.2,col=NA) +
  geom_line(size=0.3) +
  scale_y_continuous(name="Log2 spatial co-occurrence") +
  scale_x_log10(name="Distance (µm)",breaks=c(1,10),limits=c(.3,10)) +
  facet_grid(from~to) +
  annotation_logticks(sides = "b",size=0.2,
                      short = unit(0.05, "cm"),
                      mid = unit(0.1, "cm"),
                      long = unit(0.15, "cm")) +
  geom_hline(yintercept = 0,size=0.3) +
  scale_fill_manual(values = getcolor_perturbation[levels(to_plot$treatment)]) +
  scale_color_manual(values = getcolor_perturbation[levels(to_plot$treatment)]) +
  theme_bw(base_size = 7) +
  theme(legend.title = element_blank(),
        legend.key.size = unit(3,"mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0,"mm"),
        panel.grid = element_blank(),
        panel.spacing = unit(1,"mm"))
ggsave(file.path(plot_dir,"184A1_co_occurrence_Meayamycin_all.pdf"),width=12,height=12.5,units="cm")

# plot selected for main figure
to_plot <- co_occurrences_mean %>%
  # focus on Meayamycin and Unperturbed
  filter(treatment %in% c("Meayamycin (12.5h)","Unperturbed")) %>%
  filter(from!="Extra-nuclear" & to !="Extra-nuclear") %>%
  mutate(treatment=factor(treatment)) %>%
  filter(from=="Nuclear speckles" & to %in% c("Nuclear speckles","PML bodies","Nucleolus")) %>%
  mutate(to = factor(to,levels=c("Nuclear speckles","PML bodies","Nucleolus")))

# use dummy data to make symmetric resacled axes for faceting
dummy_data <- to_plot %>%
  group_by(from,to) %>%
  slice_max(abs(mean_log2_co_occurrence))

dummy_data <- bind_rows(dummy_data,mutate(dummy_data,mean_log2_co_occurrence=-mean_log2_co_occurrence))

to_plot %>%
  ggplot(aes(x=min_distance_um,y=mean_log2_co_occurrence,col=treatment,fill=treatment)) + 
  geom_ribbon(aes(ymin=mean_log2_co_occurrence_ci_lower,
                  ymax=mean_log2_co_occurrence_ci_upper),alpha=0.25,col=NA) +
  geom_line(size=0.3) +
  geom_blank(data=dummy_data) +
  scale_y_continuous(name="Log2 spatial co-occurrence\n(from Nuclear speckles)") +
  scale_x_log10(name="Distance (µm)",breaks=c(1,10),limits=c(.3,10)) +
  facet_wrap(~to,scales="free_y") +
  annotation_logticks(sides = "b",size=0.2,
                      short = unit(0.05, "cm"),
                      mid = unit(0.1, "cm"),
                      long = unit(0.15, "cm")) +
  geom_hline(yintercept = 0,size=0.3) +
  scale_fill_manual(values = getcolor_perturbation[levels(to_plot$treatment)],labels=c("Meayamycin","Unperturbed")) +
  scale_color_manual(values = getcolor_perturbation[levels(to_plot$treatment)],labels=c("Meayamycin","Unperturbed")) +
  theme_bw(base_size = 7) +
  theme(legend.title = element_blank(),
        legend.key.size = unit(3,"mm"),
        panel.grid = element_blank(),
        panel.spacing = unit(1,"mm"),
        legend.position = "bottom",
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-2,-2,-2,-2))
ggsave(file.path(plot_dir,"184A1_co_occurrence_Meayamycin_from_speckles_selected.pdf"),width=8,height=3.8,units="cm")
ggsave_cairo(file.path(plot_dir,"184A1_co_occurrence_Meayamycin_from_speckles_selected.png"),width=8,height=3.8,units="cm")


# Object statistics

boxplot_theme <- theme_bw(base_size = 7) +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        legend.direction = "vertical",
        legend.key.size = unit(2,"mm"),
        legend.box.spacing = unit(0,"mm"),
        panel.grid=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=6))

speckles_size <- object_stats %>%
  pivot_longer(cols=matches("\\|"),names_to = c("statistic","cluster_annotation"),names_sep="\\|") %>%
  filter(cluster_annotation=="Nuclear speckles" & statistic=="area_median") %>%
  filter(treatment %in% c("Meayamycin (12.5h)","Unperturbed")) %>%
  mutate(treatment=factor(treatment,levels=c("Unperturbed","Meayamycin (12.5h)"))) %>%
  mutate(value = value*pixel_size*pixel_size) 

p_speckles_size <- speckles_size %>%
  ggplot(aes(x=NA,y=value,fill=treatment)) + 
  geom_boxplot(outlier.shape = NA,size=0.25) +
  scale_y_continuous(name="Median Nuclear\nspeckle area (µm )",limits = c(0,5),breaks=scales::pretty_breaks(n=3)) +
  scale_fill_manual(values = getcolor_perturbation[levels(speckles_size$treatment)]) +
  boxplot_theme

speckles_count <- object_stats %>%
  pivot_longer(cols=matches("\\|"),names_to = c("statistic","cluster_annotation"),names_sep="\\|") %>%
  filter(cluster_annotation=="Nuclear speckles" & statistic=="count") %>%
  filter(treatment %in% c("Meayamycin (12.5h)","Unperturbed")) %>%
  mutate(treatment=factor(treatment,levels=c("Unperturbed","Meayamycin (12.5h)")))

p_speckles_count <- speckles_count %>%
  ggplot(aes(x=NA,y=value,fill=treatment)) + 
  geom_boxplot(outlier.shape = NA,size=0.25) +
  scale_y_continuous(name="Median Nuclear\nspeckle count",limits = c(0,25),breaks=scales::pretty_breaks(n=3)) +
  scale_fill_manual(values = getcolor_perturbation[levels(speckles_size$treatment)]) +
  boxplot_theme 

(p_speckles_size / p_speckles_count) +
  plot_layout(guides = "collect") &
  boxplot_theme

ggsave(file.path(plot_dir,"Meayamycin_speckles_stats.pdf"),
       width=3.1,height=4.4,units="cm")

# CX5461 ----

# plot co-occurrences
to_plot <- co_occurrences_mean %>%
  # focus on CX-5461 and Unperturbed
  filter(treatment %in% c("CX5461 (2.5h)","Unperturbed")) %>%
  filter(from!="Extra-nuclear" & to !="Extra-nuclear") %>%
  mutate(treatment=factor(treatment))

to_plot %>%
  filter(from!="Extra-nuclear" & to !="Extra-nuclear") %>%
  ggplot(aes(x=min_distance_um,y=mean_log2_co_occurrence,col=treatment,fill=treatment)) + 
  geom_ribbon(aes(ymin=mean_log2_co_occurrence_ci_lower,
                  ymax=mean_log2_co_occurrence_ci_upper),alpha=0.2,col=NA) +
  geom_line(size=0.3) +
  scale_y_continuous(name="Log2 spatial co-occurrence") +
  scale_x_log10(name="Distance (µm)",breaks=c(1,10),limits=c(.3,10)) +
  facet_grid(from~to) +
  annotation_logticks(sides = "b",size=0.2,
                      short = unit(0.05, "cm"),
                      mid = unit(0.1, "cm"),
                      long = unit(0.15, "cm")) +
  geom_hline(yintercept = 0,size=0.3) +
  scale_fill_manual(values = getcolor_perturbation[levels(to_plot$treatment)]) +
  scale_color_manual(values = getcolor_perturbation[levels(to_plot$treatment)]) +
  theme_bw(base_size = 7) +
  theme(legend.title = element_blank(),
        legend.key.size = unit(3,"mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0,"mm"),
        panel.grid = element_blank(),
        panel.spacing = unit(1,"mm"))
ggsave(file.path(plot_dir,"184A1_co_occurrence_CX5461_all.pdf"),width=12,height=12.5,units="cm")


# plot selected for main figure
to_plot <- co_occurrences_mean %>%
  # focus on CX5461 and Unperturbed
  filter(treatment %in% c("CX5461 (2.5h)","Unperturbed")) %>%
  filter(from!="Extra-nuclear" & to !="Extra-nuclear") %>%
  mutate(treatment=factor(treatment)) %>%
  filter(from=="Nucleolus" & to %in% c("Nuclear periphery","Nucleolus")) %>%
  mutate(to = factor(to,levels=c("Nucleolus","Nuclear periphery")))

# use dummy data to make symmetric resacled axes for faceting
dummy_data <- to_plot %>%
  group_by(from,to) %>%
  slice_max(abs(mean_log2_co_occurrence))

dummy_data <- bind_rows(dummy_data,mutate(dummy_data,mean_log2_co_occurrence=-mean_log2_co_occurrence))

to_plot %>%
  ggplot(aes(x=min_distance_um,y=mean_log2_co_occurrence,col=treatment,fill=treatment)) + 
  geom_ribbon(aes(ymin=mean_log2_co_occurrence_ci_lower,
                  ymax=mean_log2_co_occurrence_ci_upper),alpha=0.25,col=NA) +
  geom_line(size=0.3) +
  geom_blank(data=dummy_data) +
  scale_y_continuous(name="Log2 spatial co-occurrence\n(from Nucleolus)") +
  scale_x_log10(name="Distance (µm)",breaks=c(1,10),limits=c(.3,10)) +
  facet_wrap(~to,scales="free_y") +
  annotation_logticks(sides = "b",size=0.2,
                      short = unit(0.05, "cm"),
                      mid = unit(0.1, "cm"),
                      long = unit(0.15, "cm")) +
  geom_hline(yintercept = 0,size=0.3) +
  scale_fill_manual(values = getcolor_perturbation[levels(to_plot$treatment)],labels=c("Meayamycin","Unperturbed")) +
  scale_color_manual(values = getcolor_perturbation[levels(to_plot$treatment)],labels=c("Meayamycin","Unperturbed")) +
  theme_bw(base_size = 7) +
  theme(legend.title = element_blank(),
        legend.key.size = unit(3,"mm"),
        panel.grid = element_blank(),
        panel.spacing = unit(1,"mm"),
        legend.position = "bottom",
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-2,-2,-2,-2))
ggsave(file.path(plot_dir,"184A1_co_occurrence_CX5461_from_nucleolus_selected.pdf"),width=5.5,height=3.8,units="cm")
ggsave_cairo(file.path(plot_dir,"184A1_co_occurrence_CX5461_from_nucleolus_selected.png"),width=5.5,height=3.8,units="cm")

# Area between curves ----

area_under_curve <- function(x,y) {
  n <- length(x)
  m <- length(y)
  if (n!=m) stop("Unequal input vector lengths")
  integral <- 0
  for (i in 1:(n-1)) {
    integral <- integral + (x[i+1] - x[i]) * abs(y[i] + y[i+1]) / 2.0
  }
  return(integral)
}

area_between_curves <- co_occurrences_mean %>%
  filter(treatment %in% c("CX5461 (2.5h)","Unperturbed")) %>%
  mutate(log_min_distance_um = log(min_distance_um), log_max_distance_um = log(max_distance_um)) %>%
  group_by(from,to,treatment) %>%
  summarise(area_under_curve = area_under_curve(min_distance_um,mean_log2_co_occurrence),
            area_under_curve_log = area_under_curve(log_min_distance_um,mean_log2_co_occurrence)) %>%
  ungroup() %>%
  select(-area_under_curve) %>%
  pivot_wider(names_from = treatment,values_from = area_under_curve_log) %>%
  mutate(area_between_curves_log = abs(Unperturbed - `CX5461 (2.5h)`))

# set heatmap options
ht_opt("heatmap_row_names_gp" = gpar(fontsize=6),
       "heatmap_column_names_gp" = gpar(fontsize=6),
       "legend_title_gp" = gpar(fontsize=6,fontface="bold"),
       "legend_labels_gp" = gpar(fontsize=6),
       "legend_grid_height" = unit(2, "mm"),
       "legend_grid_width" = unit(1.5, "mm"))

hm <- area_between_curves %>%
  filter(from!="Extra-nuclear" & to !="Extra-nuclear") %>%
  select(from,to,area_between_curves_log) %>%
  pivot_wider(names_from = to,values_from = area_between_curves_log) %>%
  column_to_rownames("from") %>%
  data.matrix() %>%
  ComplexHeatmap::Heatmap(name = "Area\nbetween\ncurves",
                          col=viridisLite::magma(100),
                          column_dend_gp = gpar(lwd=0.3),
                          row_dend_gp = gpar(lwd=0.3),
                          column_dend_height = unit(0.5,"cm"),
                          row_dend_width = unit(0.5,"cm")) 
pdf(file = file.path(plot_dir,"cx5461_versus_control_area_between_curves.pdf"),
    width=6/cm(1),height=5/cm(1))
draw(hm)
dev.off()
png(file = file.path(plot_dir,"cx5461_versus_control_area_between_curves.png"),
    width=6/cm(1),height=5/cm(1),units = "in",res=300)
draw(hm)
dev.off()





