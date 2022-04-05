source_python(file.path(campa_ana$constants$SOURCE_DIR,"NascentRNA_constants.py"))

experiment_dir <- file.path(campa$constants$EXPERIMENT_DIR,experiment_name,"aggregated","full_data")
campa_res_dirs <- list.dirs(experiment_dir) %>%
  tibble(campa_res_dir=.) %>%
  filter(grepl(cell_type,campa_res_dir)) %>%
  mutate(well_name = basename(campa_res_dir)) %>%
  filter(grepl("^[A-Z]\\d{2}$",well_name)) %>%
  # check if campa_res_dir are empty
  mutate(n_files = map_int(campa_res_dir,~length(list.files(.)))) %>%
  filter(n_files > 1)

data_dirs <- list.dirs(DATA_DIR) %>%
  tibble(data_dir=.) %>%
  filter(grepl(cell_type,data_dir)) %>%
  mutate(well_name = basename(data_dir)) %>%
  filter(grepl("^[A-Z]\\d{2}$",well_name)) %>%
  filter(!grepl("ilastik",data_dir)) %>%
  # check if data_dir are empty
  mutate(n_files = map_int(data_dir,~length(list.files(.)))) %>%
  filter(n_files > 1)

ilastik_data_dirs <- list.dirs(DATA_DIR) %>%
  tibble(ilastik_data_dir=.) %>%
  filter(grepl(cell_type,ilastik_data_dir)) %>%
  mutate(well_name = basename(ilastik_data_dir)) %>%
  filter(grepl("^[A-Z]\\d{2}$",well_name)) %>%
  filter(grepl("ilastik",ilastik_data_dir)) %>%
  # check if ilastik_data_dir are empty
  mutate(n_files = map_int(ilastik_data_dir,~length(list.files(.)))) %>%
  filter(n_files > 1)

