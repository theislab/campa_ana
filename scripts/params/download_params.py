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
        'well_names': ['I11'], #['I09', 'I11', 'J10', 'J12'],
        # output directory name (appended to DATA_DIR)
        'out_dir': '184A1_unperturbed',
        # Cells or Nuclei
        'object': 'Nuclei',
        # channel names that should be downloaded
        'channels': ['00_EU'], #download_channels['default'],
        # mean object size to estimate allocated memory
        'size': 50000
    },
    'AZD4573': {
        'well_names': ['I13', 'I17', 'J14', 'J18', 'J21'],
        # output directory name (appended to DATA_DIR)
        'out_dir': '184A1_AZD4573',
        # Cells or Nuclei
        'object': 'Nuclei',
        # channel names that should be downloaded
        'channels': download_channels['default'],
        # mean object size to estimate allocated memory
        'size': 50000
    } # TODO add remaining params
}