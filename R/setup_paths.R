source_python(file.path(campa_ana$constants$SOURCE_DIR,"NascentRNA_constants.py"))

experiment_dir <- file.path(campa$constants$EXPERIMENT_DIR,experiment_name,"aggregated","full_data")
campa_res_dirs <- list.dirs(experiment_dir) %>%
  tibble(campa_res_dir=.) %>%
  filter(grepl(cell_type,campa_res_dir)) %>%
  mutate(well_name = basename(campa_res_dir)) %>%
  filter(grepl("^[A-Z]\\d{2}$",well_name)) 

data_dirs <- list.dirs(DATA_DIR) %>%
  tibble(data_dir=.) %>%
  filter(grepl(cell_type,data_dir)) %>%
  mutate(well_name = basename(data_dir)) %>%
  filter(grepl("^[A-Z]\\d{2}$",well_name)) 

