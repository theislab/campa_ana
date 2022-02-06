import pandas as pd
import os
import argparse
import mcu
import mcu.mpp
from campa_ana.tm_data import get_tm_client
from campa.constants import get_data_config
from campa.utils import load_config

def download_mpp(args):
    params = load_config(args.params)
    tm_credentials = params.tm_credentials
    experiment_name = params.experiment_name
    download_params = params.download_params[args.key]
    tm_client = get_tm_client(tm_credentials, experiment_name)
    DATA_CONFIG = get_data_config(params.data_config)
    
    for well in download_params['well_names']:
        cur_dir = os.path.join(DATA_CONFIG.DATA_DIR, download_params['out_dir'], well)
        if not os.path.exists(cur_dir):
            os.makedirs(cur_dir)
        # get metadata
        metadata_nuclei = tm_client.download_object_metadata('Nuclei', well_name=well)
        metadata_cells = tm_client.download_object_metadata('Cells', well_name=well)
        # merge metadata to get cell mapoject ids
        metadata = metadata_nuclei.join(metadata_cells, how='left', rsuffix='_cell')
        # make sure that labels match
        for col in ['label', 'plate_name', 'well_name', 'well_pos_x', 'well_pos_y', 'tpoint', 'zplane']:
            assert (metadata[col] == metadata[col+'_cell']).all(), col
            metadata.drop(columns=col+'_cell')

        selected_metadata = metadata[metadata.is_border==0]
        # save metadata
        selected_metadata.to_csv(os.path.join(cur_dir, 'metadata.csv'))
        
        mpp_args = argparse.Namespace(metadata_file=os.path.join(cur_dir, 'metadata.csv'), 
                                  credentials_file=tm_credentials,
                                  output_directory=cur_dir,
                                  object_type=download_params['object'],
                                  experiment_name=experiment_name,
                                  channel_names=download_params['channels'],
                                  mean_object_size=download_params['size'])
        #mcu.mpp.main(mpp_args)

def parse_arguments():
    parser = argparse.ArgumentParser(
        description=('Download MPPs from different wells')
    )
    parser.add_argument('--params', help='download_params file', default='download_mpp_params.py')
    parser.add_argument('--key', help='key in params dict identifying specific params list that should be downloaded', default='unperturbed')
    parser.add_argument('--verbose', '-v', action='count', default=0)
    return(parser.parse_args())


if __name__ == "__main__":
    args = parse_arguments()
    mcu.mpp.setup_logger(args)
    download_mpp(args)