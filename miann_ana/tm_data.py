from tmclient import TmClient
import pandas as pd
import os

def get_tm_client(credentials_file, experiment_name='201908-NascentRNA-4i'):
    """Returns TmClient instance from tm_credentials file (csv)"""
    tm_credentials = pd.read_csv(credentials_file)
    tm_client = TmClient(
        host=tm_credentials.host[0],
        port=tm_credentials.port[0],
        username=tm_credentials.username[0],
        password=tm_credentials.password[0],
        experiment_name=experiment_name
    )
    return tm_client 

  
def download_channels_metadata(tm_client, output_file):
    """create csv file contraining metadata about the channels
    """
    # get data from TM
    channels = tm_client.get_channels()
    channels_data = {
        'id': [ch['id'] for ch in channels],
        'name': [ch['name'] for ch in channels],
        'max_intensity': [ch['layers'][0]['max_intensity'] for ch in channels],
        'min_intensity': [ch['layers'][0]['min_intensity'] for ch in channels]
    }
    # save to file
    pd.DataFrame(channels_data).to_csv(output_file)
    
def download_wells_metadata(tm_client, output_file, data_dir):
    # get data from TM
    wells = pd.DataFrame(tm_client.get_wells())

    # read perturbation and cell type info from Scotts well_names.csv file
    well_info = pd.read_csv(os.path.join(data_dir, 'metadata_raw/well_names.csv'))
    # merge info
    wells = wells.merge(well_info, left_on='name', right_on='well_name')
    wells = wells.rename(columns={'cell_line': 'cell_type', 'compound': 'perturbation'})
    wells.perturbation.fillna('normal', inplace=True)
    # drop "name" column (identical to "well_name")
    wells = wells.drop(columns=['name'])
    # add perturbation_duration column
    duration_str = '-' + wells['duration'].fillna(0).astype(int).astype(str)
    duration_str[wells['duration'].isna()] = ''
    wells['perturbation_duration'] = wells['perturbation'] + duration_str

    # save to file
    wells.to_csv(output_file)