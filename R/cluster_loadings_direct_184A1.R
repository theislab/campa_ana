# clear the workspace
rm(list=ls())

# load required packages
library(tidyverse)
library(patchwork)
library(ComplexHeatmap)

# setup script-specific parameters
cell_type <- "184A1"
experiment_name <- "VAE_all/MPPleiden"

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
plot_dir <- file.path(campa_ana$constants$SOURCE_DIR,"figures","cluster_loadings")
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
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

resolutions <- c("0.2","0.4","0.6","0.8","1.2","1.6","2.0")

# load campa results (per cell)
campa_leiden_res <- lapply(X = resolutions, function (x) {
  selected_wells %>%
    mutate(model_res = map(campa_res_dir,~read_csv(file.path(.,"export",paste0("intensity_features_clustering_res",x,".csv")),show_col_types = FALSE) %>% select(-1))) %>%
    select(-well_name) %>%
    mutate(leiden_resolution=as.numeric(x)) %>%
    unnest(model_res)
}
)

campa_leiden_res <- bind_rows(campa_leiden_res)

column_labels <- structure(rename_channels_for_plotting_dict$new, names = rename_channels_for_plotting_dict$old)

total_pixel_sizes <- campa_leiden_res %>%
  group_by(cluster,leiden_resolution,treatment) %>%
  summarise(size=sum(size)) %>%
  ungroup()

relative_pixel_sizes <- total_pixel_sizes %>% 
  filter(cluster!="all") %>%
  left_join(filter(total_pixel_sizes,cluster=="all") %>% 
              rename(all=size) %>%
              select(-cluster)) %>%
  mutate(relative_size=size/all) %>%
  select(-size,-all)

# Make plots ----
# PLOT - show mapping of leiden clusters to compartments using loading matrix
for (res in as.numeric(resolutions)) {
  
  by_leiden <- campa_leiden_res %>%
    filter(leiden_resolution==res) %>%
    # convery to sum intensity to average across all pixels rather than averaging per-cell values
    mutate(across(matches("^\\d{2}_"),~.*size)) %>%
    group_by(cluster) %>%
    summarise(across(matches("^\\d{2}_|size"),sum)) %>%
    mutate(across(matches("^\\d{2}_"),~./size)) %>%
    select(-size) %>%
    filter(cluster!="all") %>%
    column_to_rownames("cluster") %>%
    data.matrix() %>%
    scale() 
  
  order1 <- seriation::seriate(dist(by_leiden),method="OLO_ward")
  order2 <- seriation::seriate(dist(t(by_leiden)),method="OLO_ward")
  
  leiden_hm <- ComplexHeatmap::Heatmap(
    by_leiden,
    width = ncol(by_leiden)*unit(3, "mm"),
    height = nrow(by_leiden)*unit(3, "mm"),
    name="Mean\nIntensity\n(z-score)",
    col = cols_pm(3),
    column_labels = column_labels[colnames(by_leiden)],
    cluster_columns = as.dendrogram(order2[[1]]),
    cluster_rows = as.dendrogram(order1[[1]]),
    row_dend_side = "left",
    row_names_side = "right",
    #right_annotation = ha,
    border = "black",
    row_title_gp = gpar(fontsize=8),
    row_title_rot = 0,
    row_names_gp = gpar(fontsize=8),
    column_names_gp = gpar(fontsize=8),
    heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "bold"),
                                labels_gp = gpar(fontsize = 8))
  )
  draw(leiden_hm)
  pdf(file = file.path(plot_dir,paste0("mpp_cluster_loadings_resolution_",res,".pdf")),width=14/cm(1),height=14/cm(1))
  draw(leiden_hm,heatmap_legend_side = "left")
  dev.off()
  png(file = file.path(plot_dir,paste0("mpp_cluster_loadings_resolution_",res,".png")),width=14/cm(1),height=14/cm(1),units="in",res=600)
  draw(leiden_hm,heatmap_legend_side = "left")
  dev.off()
  
  # Now compute changes in cluster size across conditions
  log2_fc <- relative_pixel_sizes %>% 
    filter(treatment!="Control" & leiden_resolution==res) %>%
    left_join(filter(relative_pixel_sizes,treatment=="Control") %>% 
                rename(Control=relative_size) %>%
                select(-treatment)) %>%
    mutate(log2FC_relative_size=log2(relative_size/Control)) %>%
    select(-Control) 
  
  log2_fc_mat <- log2_fc %>%
    select(cluster,treatment,log2FC_relative_size) %>%
    pivot_wider(names_from = treatment, values_from = log2FC_relative_size) %>%
    column_to_rownames("cluster") %>%
    data.matrix() 
  
  cluster_order <- seriation::get_order(order1)
  treatment_order <- seriation::get_order(seriation::seriate(dist(t(log2_fc_mat)),method="OLO_ward"))
  
  rownames(log2_fc_mat)[cluster_order]
  colnames(log2_fc_mat)[treatment_order]
  
  log2_fc_plot <- log2_fc %>%
    mutate(cluster=factor(cluster,levels=rev(rownames(log2_fc_mat)[cluster_order]))) %>%
    mutate(treatment=factor(treatment,levels=colnames(log2_fc_mat)[treatment_order])) %>%
    ggplot(aes(x=cluster,y=treatment)) + 
    geom_point(aes(size=relative_size,fill=log2FC_relative_size),stroke=0.2,pch=21) +
    scale_fill_gradient2(name=expression(log[2](`fold-change`)),
                         breaks=c(-5,0,5),
                         low = "#742a7e",mid="white",high = "#1a572b",midpoint = 0,
                         limits=c(-5,5),oob=scales::squish,
                         guide = guide_colourbar(order=1)) +
    xlab("Cluster") +
    scale_radius(name = "Fraction of pixels",
                 range = c(1.5,2.8),
                 limits = c(0,NA),
                 guide = guide_legend(order=2)) +
    theme_bw(base_size = 7) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle=90,hjust=1,vjust=0.5),
          axis.text.y = element_text(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none",
          aspect.ratio = nrow(log2_fc_mat)/ncol(log2_fc_mat)) + coord_flip()
  log2_fc_plot 
  
  h <- 2.35 + 0.225*nrow(log2_fc_mat)
  ggsave(plot = log2_fc_plot,filename = file.path(plot_dir,paste0("direct_clustering_log2FC_dot_plot_res_",res,".pdf")),width=3.5,height=h,units="cm")
  ggsave_cairo(plot = log2_fc_plot,filename = file.path(plot_dir,paste0("direct_clustering_log2FC_dot_plot_res_",res,".png")),width=3.5,height=h,units="cm",dpi=600)
  
  p <- log2_fc_plot + theme(
    legend.key.size = unit(0.15,"cm"),
    legend.title = element_text(size=6),
    legend.box = "vertical",
    legend.position="right",
    legend.box.margin = margin(0,0,0,0,unit="cm")) 
  
  legend_plot <- ggpubr::get_legend(p) %>%
    ggpubr::as_ggplot()
  legend_plot
  ggsave(plot = legend_plot,filename = file.path(plot_dir,paste0("direct_clustering_log2FC_dot_plot_res_",res,"_legend.pdf")),width=3.5,height=3.5,units="cm")
  ggsave_cairo(plot = log2_fc_plot,filename = file.path(plot_dir,paste0("direct_clustering_log2FC_dot_plot_res_",res,"_legend.png")),width=3.5,height=3.5,units="cm",dpi=600)
  
}

# add colorbar to resolution 1.2 plot
res1_2_colors <- c("#1f77b4","#aec7e8","#e377c2","#f7b6d2","#7f7f7f","#c7c7c7","#bcbd22","#dbdb8d","#17becf",
"#9edae5","#ff7f0e","#ffbb78","#2ca02c","#98df8a","#d62728","#ff9896","#9467bd","#c5b0d5",
"#8c564b","#c49c94","#ffffff")
res1_2_colors <- setNames(res1_2_colors, as.character(0:20))

res <- "1.2"

by_leiden <- campa_leiden_res %>%
  filter(leiden_resolution==res) %>%
  # convery to sum intensity to average across all pixels rather than averaging per-cell values
  mutate(across(matches("^\\d{2}_"),~.*size)) %>%
  group_by(cluster) %>%
  summarise(across(matches("^\\d{2}_|size"),sum)) %>%
  mutate(across(matches("^\\d{2}_"),~./size)) %>%
  select(-size) %>%
  filter(cluster!="all") %>%
  column_to_rownames("cluster") %>%
  data.matrix() %>%
  scale() 

order1 <- seriation::seriate(dist(by_leiden),method="OLO_ward")
order2 <- seriation::seriate(dist(t(by_leiden)),method="OLO_ward")

ha = HeatmapAnnotation(
  #foo = as.character(seriation::get_order(order1)),
  foo = sort(as.character(0:19)),
  col = list(foo = res1_2_colors),
  na_col = "black",
  which = "row",
  show_annotation_name = FALSE,
  show_legend = FALSE
)

leiden_hm <- ComplexHeatmap::Heatmap(
  by_leiden,
  width = ncol(by_leiden)*unit(3, "mm"),
  height = nrow(by_leiden)*unit(3, "mm"),
  name="Mean\nIntensity\n(z-score)",
  col = cols_pm(3),
  column_labels = column_labels[colnames(by_leiden)],
  cluster_columns = as.dendrogram(order2[[1]]),
  cluster_rows = as.dendrogram(order1[[1]]),
  row_dend_side = "left",
  row_names_side = "right",
  right_annotation = ha,
  border = "black",
  row_title_gp = gpar(fontsize=8),
  row_title_rot = 0,
  row_names_gp = gpar(fontsize=8),
  column_names_gp = gpar(fontsize=8),
  heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "bold"),
                              labels_gp = gpar(fontsize = 8))
)
draw(leiden_hm)
pdf(file = file.path(plot_dir,paste0("mpp_cluster_loadings_resolution_",res,"_with_annotation.pdf")),width=14/cm(1),height=14/cm(1))
draw(leiden_hm,heatmap_legend_side = "left")
dev.off()
png(file = file.path(plot_dir,paste0("mpp_cluster_loadings_resolution_",res,"_with_annotation.png")),width=14/cm(1),height=14/cm(1),units="in",res=600)
draw(leiden_hm,heatmap_legend_side = "left")
dev.off()

