from campa.constants import campa_config
from campa.tl import Experiment, Predictor, Cluster
from campa_ana.constants import SOURCE_DIR
from campa.data import MPPData
import os
from campa.tl._cluster import annotate_clustering
import numpy as np
import pandas as pd
import tensorflow as tf
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import pickle

# Integrated gradients function
# Code adapted from here: https://www.tensorflow.org/tutorials/interpretability/integrated_gradients

def interpolate_mpps(baseline,mpp,alphas):
    alphas_x = alphas[:, tf.newaxis, tf.newaxis, tf.newaxis]
    baseline_x = tf.expand_dims(baseline, axis=0)
    input_x = tf.expand_dims(mpp, axis=0)
    delta = input_x - baseline_x
    mpps = baseline_x +  alphas_x * delta
    return mpps

def compute_gradients(pred, mpp, condition, channel_idx=0, use_model='model'):
    #latent = pred.est.model.encoder((mpps, tf.tile(condition, [len(mpps),1])))
    model_input = (mpp, condition)
    #print(mpp.shape, condition.shape)
    with tf.GradientTape() as tape:
        tape.watch(model_input)
        if use_model == 'model':
            res = pred.est.model.model(model_input)[0][:,channel_idx]
        elif use_model == 'decoder':
            res = pred.est.model.decoder(model_input)[:, channel_idx]
    return tape.gradient(res, model_input)

def integral_approximation(gradients):
    # gradients is tuple
    integrated_gradients = []
    for grad in gradients:
        # riemann_trapezoidal
        grads = (grad[:-1] + grad[1:]) / tf.constant(2.0, dtype=grad.dtype)
        integrated_gradients.append(tf.math.reduce_mean(grads, axis=0))
    return integrated_gradients

def integrated_gradients(pred, mpp, condition, m_steps=50, num_channels=34, baseline_mpp=None, use_model='model'):
    if baseline_mpp is None:
        baseline_mpp = tf.zeros_like(mpp)
    #print(mpp.shape, baseline_mpp.shape)
    baseline_condition = tf.zeros_like(condition)
    # Generate alphas
    alphas = tf.linspace(start=0.0, stop=1.0, num=m_steps+1)
    # Get Gradients
    interpolated_mpps = interpolate_mpps(baseline_mpp, mpp, tf.cast(alphas, mpp.dtype))
    if use_model == 'decoder':
        interpolated_mpps = interpolated_mpps[:,0,0,:]
    interpolated_conditions = interpolate_mpps(baseline_condition, condition, tf.cast(alphas, dtype=condition.dtype))[:,0,0,:]
    print('interp mpps', interpolated_mpps.shape, 'interp conditions', interpolated_conditions.shape)
    
    integrated_gradients = []
    for channel_idx in range(num_channels):
        gradients = compute_gradients(pred, interpolated_mpps, interpolated_conditions, channel_idx=channel_idx, use_model=use_model)
        # Integral approximation through averaging gradients.
        avg_gradients = integral_approximation(gradients=gradients)

        # Scale integrated gradients with respect to input.
        integrated_gradients.append([(mpp - baseline_mpp) * avg_gradients[0], (condition - baseline_condition) * avg_gradients[1]])
        #integrated_gradients = (mpp - baseline) * avg_gradients

    return integrated_gradients

def gradXinput(pred, mpps, conditions, channel_idx=0):
    latents = pred.est.model.encoder((mpps, conditions))
    gradients = compute_gradients(pred, latents, conditions, channel_idx=channel_idx)
    return [gradients[0]*latents, gradients[1]*conditions]


if __name__ == "__main__":
    # switch to determine to what gradients are calculated. 
    # 'decoder' means grads are calculated wrt latent + condition
    # 'model' means grads are calculated wrt input + condition
    use_model = 'decoder' # model
    save_prefix = 'saliency_output_latent'  # saliency_output_v2
    fig_dir = Path(SOURCE_DIR)/'figures'/'fig1_suppl'
    os.makedirs(str(fig_dir), exist_ok=True)

    # load data
    exp = Experiment.from_dir('VAE_all/CondVAE_pert-CC')
    pred = Predictor(exp)
    mpp_data = MPPData.from_data_dir('aggregated/sub-0.001', data_config='NascentRNA', base_dir=exp.full_path, keys=['latent', 'conditions', 'clustering_res0.5'])
    cluster_annotation = pd.read_csv(os.path.join(exp.full_path, 'aggregated/sub-0.001/clustering_res0.5_annotation.csv'), index_col=0)
    # add annotation of clustering to mpp_data
    mpp_data._data['annotation'] = annotate_clustering(list(map(int, mpp_data.data('clustering_res0.5'))), annotation=cluster_annotation, cluster_name='clustering_res0.5', annotation_col='annotation')

    # calculate baseline: mean nucleoplasm intensity of unperturbed G1 cells
    obj_mask = mpp_data.metadata['perturbation_duration'].isin(['normal', 'DMSO-720', 'DMSO-120']) & (mpp_data.metadata['cell_cycle'] == 'G1')
    mpp_mask = mpp_data._get_per_mpp_value(obj_mask)
    

    nucleoplasm_mask = mpp_data._data['annotation'] == 'Nucleoplasm'
    if use_model == 'model':
        background = mpp_data.mpp[nucleoplasm_mask & mpp_mask].mean(axis=0)
        baseline=tf.convert_to_tensor(background)
    elif use_model == 'decoder':
        # need to get latent background
        print('getting baseline')
        mpp_bg = tf.convert_to_tensor(mpp_data.mpp[nucleoplasm_mask & mpp_mask])
        cond_bg = tf.convert_to_tensor(mpp_data.conditions[nucleoplasm_mask & mpp_mask])
        baseline = tf.math.reduce_mean(pred.est.model.encoder((mpp_bg, cond_bg)), axis=0)

    # calculate gradients for all samples
    perturbations = [
        ['normal', 'DMSO-720', 'DMSO-120'],
        ['AZD4573-30'], 
        ['TSA-30'], 
        ['CX5461-120'], 
        ['AZD4573-120'], 
        ['Triptolide-120'], 
        ['Meayamycin-720']
        ]

    for perts in perturbations:
        print(f'calculating saliency for perts {perts}')
        mpp_mask = mpp_data._get_per_mpp_value(mpp_data.metadata['perturbation_duration'].isin(perts))

        grads_mpp = []
        grads_condition = []
        for i in range(len(mpp_data.mpp[mpp_mask])):
            mpp = tf.convert_to_tensor(mpp_data.mpp[mpp_mask][i])
            cond = tf.expand_dims(tf.convert_to_tensor(mpp_data.conditions[mpp_mask][i]), axis=0)
            if use_model == 'decoder':
                # need to calc latent first
                print('before encoder', mpp.shape, cond.shape)
                mpp = pred.est.model.encoder((mpp[tf.newaxis], cond))[0]
                print('after encoder', mpp.shape)
                print('baseline shape', baseline.shape)
            res = integrated_gradients(pred, mpp, cond, num_channels=len(mpp_data.channels), baseline_mpp=baseline, m_steps=20, use_model=use_model)
            cur_grads_mpp = []
            cur_grads_condition = []
            for cur_res in res:
                cur_grads_mpp.append(cur_res[0].numpy())
                cur_grads_condition.append(cur_res[1].numpy())
            grads_mpp.append(cur_grads_mpp)
            grads_condition.append(cur_grads_condition)
            
        grads_mpp = np.array(grads_mpp)
        grads_condition = np.array(grads_condition)
        if use_model == 'model':
            grads_mpp = grads_mpp.transpose((1,0,2,3,4))
            grads_condition = grads_condition.transpose((1,0,2,3))

        # results shape: output channel, num data, ..., input_channel 
        pickle.dump((grads_mpp, grads_condition), open(fig_dir / f'{save_prefix}_{"_".join(perts)}.pkl', 'wb'))