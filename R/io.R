ggsave_cairo <- function(filename,plot = last_plot(),width=NA,height=NA,units="cm",dpi=300) {
  require(Cairo)
  CairoPNG(filename=filename,height=height,width=width,units=units,res=dpi)
  print(plot)
  dev.off()
}

load_channels <- function(data_dir) {
  channels <- read_csv(file.path(data_dir,"channels.csv"),col_names = F,show_col_types = F)
  names(channels) <- c("index","channel_name")
  return(channels$channel_name)
}

load_obj_ids <- function(dirname) {
  mapobject_ids <- np$load(file.path(dirname,"obj_ids.npy"))
  return(tibble(mapobject_id=mapobject_ids))
}

load_input_data <- function(dirname,channels) {
  mapobject_ids <- np$load(file.path(dirname,"obj_ids.npy"))
  labels <- np$load(file.path(dirname,"labels.npy"))
  x <- np$load(file.path(dirname,"x.npy"))
  y <- np$load(file.path(dirname,"y.npy"))
  
  mpp <- np$load(file.path(dirname,"mpp.npy"))
  colnames(mpp) <- channels
  combined <- tibble(mapobject_id=mapobject_ids,
                     x=x,
                     y=y) %>%
    bind_cols(as_tibble(mpp))
  return(combined)
}

load_ilastik_data <- function(dirname,channels) {
  mapobject_ids <- np$load(file.path(dirname,"mapobject_ids.npy"))
  labels <- np$load(file.path(dirname,"labels.npy"))
  x <- np$load(file.path(dirname,"x.npy"))
  y <- np$load(file.path(dirname,"y.npy"))
  
  mpp <- np$load(file.path(dirname,"mpp.npy"))
  colnames(mpp) <- channels
  combined <- tibble(mapobject_id=mapobject_ids,
                     x=x,
                     y=y) %>%
    bind_cols(as_tibble(mpp))
  return(combined)
}

load_pixel_clustering <- function(dirname,filename) {
  mapobject_ids <- np$load(file.path(dirname,"obj_ids.npy"))
  cluster_ids <- np$load(file.path(dirname,filename))
  x <- np$load(file.path(dirname,"x.npy"))
  y <- np$load(file.path(dirname,"y.npy"))
  
  combined <- tibble(mapobject_id=mapobject_ids,
                     x=x,
                     y=y,
                     cluster_id=cluster_ids) 
  return(combined)
}

load_latent_space <- function(dirname) {
  mapobject_ids <- np$load(file.path(dirname,"obj_ids.npy"))
  latent_space <- np$load(file.path(dirname,"latent.npy"))
  
  combined <- tibble(mapobject_id=mapobject_ids,
                     latent_space=latent_space) 
  return(combined)
}
