# --- download parameters for downloading all used NascentRNA data from tissuemaps ---
# name of data config for determining where data should be saved
data_config = "NascentRNA"
# path to tm credentials file
tm_credentials = "/Users/hannah.spitzer/projects/pelkmans/tm_credentials"
# name of the experiment in TissueMaps
experiment_name = "201908-NascentRNA-4i"

download_channels = {
    # new list of channels (2021-04)
    'default': ['00_EU', '01_CDK9_pT186', '01_PABPC1', '02_CDK7', '03_CDK9', '03_RPS6', '05_GTF2B', '05_Sm', '07_POLR2A', '07_SETD1A', '08_H3K4me3', '09_CCNT1', '09_SRRM2', '10_H3K27ac', '10_POL2RA_pS2', '11_KPNA2_MAX', '11_PML', '12_RB1_pS807_S811', '12_YAP1', '13_PABPN1', '13_POL2RA_pS5', '14_PCNA', '15_SON', '15_U2SNRNPB', '16_H3', '17_HDAC3', '17_SRSF2', '18_NONO', '19_KPNA1_MAX', '20_ALYREF', '20_SP100', '21_COIL', '21_NCL', '00_DAPI', '07_H2B'], 
    # channels for SFB2 experiments (2021-04)
    'SBF2': ['00_EU', '01_CDK9_pT186', '01_PABPC1', '02_CDK7', '02_RPS6_pS235_S236', '03_CDK9', '03_RPS6', '04_MAPK1_pT202_T204', '05_GTF2B', '05_Sm', '06_CALR', '06_CTNNB1', '07_POLR2A', '07_SETD1A', '08_H3K4me3', '08_PXN', '09_CCNT1', '09_SRRM2', '10_H3K27ac', '10_POL2RA_pS2', '11_KPNA2_MAX', '11_PML', '12_RB1_pS807_S811', '12_YAP1', '13_PABPN1', '13_POL2RA_pS5', '14_HSPD1', '14_PCNA', '15_SON', '15_U2SNRNPB', '16_GOLGA2', '16_H3', '17_HDAC3', '17_SRSF2', '18_NONO', '19_KPNA1_MAX', '19_TUBA1A', '20_ALYREF', '20_SP100', '21_COIL', '21_NCL', '22_DDX6', '00_DAPI', '07_H2B'] 
}
download_params = {
    'unperturbed': {
        'well_names': ['I09', 'I11', 'J10', 'J12'],
        # output directory name (appended to DATA_DIR)
        'out_dir': '184A1_unperturbed',
        # Cells or Nuclei
        'object': 'Nuclei',
        # channel names that should be downloaded
        'channels': download_channels['default'],
        # mean object size to estimate allocated memory
        'size': 50000
    },
    'AZD4573': {
        'well_names': ['I13', 'I17', 'J14', 'J18', 'J21'],
        'out_dir': '184A1_AZD4573',
        'object': 'Nuclei',
        'channels': download_channels['default'],
        'size': 50000
    },
    'CX5461': {
        'well_names': ['I18', 'J09', 'J22'],
        'out_dir': '184A1_CX5461',
        'object': 'Nuclei',
        'channels': download_channels['default'],
        'size': 50000
    },
    'DMSO': {
        'well_names': ['I14', 'J16'],
        'out_dir': '184A1_DMSO',
        'object': 'Nuclei',
        'channels': download_channels['default'],
        'size': 50000
    },
    'TSA': {
        'well_names': ['I16', 'J13', 'J20'],
        'out_dir': '184A1_TSA',
        'object': 'Nuclei',
        'channels': download_channels['default'],
        'size': 50000
    },
    'triptolide': {
        'well_names': ['I10', 'J15'],
        'out_dir': '184A1_triptolide',
        'object': 'Nuclei',
        'channels': download_channels['default'],
        'size': 50000
    },
    'meayamycin': {
        'well_names': ['I12', 'I20'],
        'out_dir': '184A1_meayamycin',
        'object': 'Nuclei',
        'channels': download_channels['default'],
        'size': 50000
    },
    # EU10 wells with 10 mins EU for conditioning on TR
    'EU10': {
        'well_names': ['H07', 'I08', 'J07'],
        'out_dir': '184A1_EU10',
        'object': 'Nuclei',
        'channels': download_channels['default'],
        'size': 50000
    },
    # siRNA wells + controls for training on entire cells
    'SBF2': {
        'well_names': ['K18', 'L19', 'M18'],
        'out_dir': '184A1_SBF2',
        'object': 'Cells',
        'channels': download_channels['SBF2'],
        'size': 500000  # TODO check on pelkmanslab server what size was used
    },
    'scrambled': {
        'well_names': ['K19', 'L18', 'M19'],
        'out_dir': '184A1_scrambled',
        'object': 'Cells',
        'channels': download_channels['SBF2'],
        'size': 500000
    },
}