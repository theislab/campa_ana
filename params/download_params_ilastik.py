# --- download parameters for downloading all used NascentRNA data from tissuemaps ---
# name of data config for determining where data should be saved
data_config = "NascentRNA"
# path to tm credentials file
#tm_credentials = "/Users/hannah.spitzer/projects/pelkmans/tm_credentials"
tm_credentials = "/data/homes/hannah/NascentRNA/tm_credentials"
# name of the experiment in TissueMaps
experiment_name = "201908-NascentRNA-4i"

download_channels = {
    # list of ilastik channels
    'default': ['09_SRRM2_ILASTIK', '11_PML_ILASTIK', '15_SON_ILASTIK', '20_SP100_ILASTIK'],
   }
download_params = {
    'unperturbed': {
        'well_names': ['I09', 'I11', 'J10', 'J12'],
        # output directory name (appended to DATA_DIR)
        'out_dir': 'ilastik/184A1_unperturbed',
        # Cells or Nuclei
        'object': 'Nuclei',
        # channel names that should be downloaded
        'channels': download_channels['default'],
        # mean object size to estimate allocated memory
        'size': 50000,
        # illumination correction? (default is True)
        'correct': False
    },
    'AZD4573': {
        'well_names': ['I13', 'I17', 'J14', 'J18', 'J21'],
        'out_dir': 'ilastik/184A1_AZD4573',
        'object': 'Nuclei',
        'channels': download_channels['default'],
        'size': 50000,
        'correct': False,
    },
    'CX5461': {
        'well_names': ['I18', 'J09', 'J22'],
        'out_dir': 'ilastik/184A1_CX5461',
        'object': 'Nuclei',
        'channels': download_channels['default'],
        'size': 50000,
        'correct': False
    },
    'DMSO': {
        'well_names': ['I14', 'J16'],
        'out_dir': 'ilastik/184A1_DMSO',
        'object': 'Nuclei',
        'channels': download_channels['default'],
        'size': 50000,
        'correct': False
    },
    'TSA': {
        'well_names': ['I16', 'J13', 'J20'],
        'out_dir': 'ilastik/184A1_TSA',
        'object': 'Nuclei',
        'channels': download_channels['default'],
        'size': 50000,
        'correct': False
    },
    'triptolide': {
        'well_names': ['I10', 'J15'],
        'out_dir': 'ilastik/184A1_triptolide',
        'object': 'Nuclei',
        'channels': download_channels['default'],
        'size': 50000,
        'correct': False
    },
    'meayamycin': {
        'well_names': ['I12', 'I20'],
        'out_dir': 'ilastik/184A1_meayamycin',
        'object': 'Nuclei',
        'channels': download_channels['default'],
        'size': 50000,
        'correct': False
    },
    # EU10 wells with 10 mins EU for conditioning on TR
    'EU10': {
        'well_names': ['H07', 'I08', 'J07'],
        'out_dir': 'ilastik/184A1_EU10',
        'object': 'Nuclei',
        'channels': download_channels['default'],
        'size': 50000,
        'correct': False
    },
}