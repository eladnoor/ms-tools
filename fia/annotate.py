#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 23:54:03 2015

@author: eladn
"""

import numpy as np
import sys
import os
import argparse
import pandas as pd
from tqdm import tqdm

base_path = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(base_path)
from openbis import download_data_profiles, get_sample_names

# script parameters (as determined by Mattia)
MIN_PEAK_SIZE = 5000
MAX_MZ_DIFFERENCE = 0.003
REF_MASS_RANGE = (50, 1000)
REFERENCE_MASS_FNAME = os.path.join(base_path, 'EMDTB.csv')
if not os.path.exists(REFERENCE_MASS_FNAME):
    raise Exception('Cannot locate the CSV file containing reference masses: '
                    + REFERENCE_MASS_FNAME)

def findpeaks(a):
    """
    Returns:
        list of indices of the local maxima in a 1D array 'a'. A local peak is 
        larger than its two neighboring samples. Endpoints are excluded.
        If a peak is flat, the function returns only the point with the lowest index.
    """
    tmp = np.diff(np.sign(np.diff(a.flat)))
    return np.where(tmp == -2)[0] + 1

#%%
parser = argparse.ArgumentParser(description='Download FIA raw data from openBIS')
parser.add_argument('exp_code', type=str,
                    help='the openBIS experiment ID')
parser.add_argument('-o', dest='output_path', type=str, default=None,
                    help='the path where all output files will be written')
args = parser.parse_args()

#%%
sample_df = pd.Series(get_sample_names(args.exp_code))
sample_df.index.name = 'sample.code'
sample_df.name = 'sample.name'
sample_df = sample_df.to_frame().sort_index()

dataProfiles = download_data_profiles(args.exp_code)
dsSampleCodes = sorted(dataProfiles.keys())
n_samples = len(dsSampleCodes)

# allPeaks is a list of matrices (one per sample) of the ion mz and intensity
# only for the peaks (i.e. local maxima)
allPeaks = {}

#%% identify peaks (local maxima)
for s in tqdm(dsSampleCodes, desc='Identifying centroids'):

    # find all the values that are local maxima and pass the threshold
    idxs = findpeaks(dataProfiles[s][:, 1])
    idxs = list(filter(lambda j : dataProfiles[s][j, 1] >= MIN_PEAK_SIZE, idxs))
    allPeaks[s] = dataProfiles[s][idxs, :]

#%% Use the reference table to associate peaks to compounds, by minimum mass distance
reference_df = pd.read_csv(REFERENCE_MASS_FNAME)
reference_df.index.name = 'ion.code'

# subtract the mass of H+ (i.e. look for deprotonated masses)
proton_mass = reference_df.loc[0, 'mass']

# keep only ions in the relevant range for FIA
compound_df = reference_df[(REF_MASS_RANGE[0] < reference_df['mass']) &
                           (REF_MASS_RANGE[1] > reference_df['mass'])]

# peak_masses[i, j] will contain the exact mass of the peak which is closest
# to reference 'j' in sample 'i'. If there is no peak which is close enough
# (i.e. in the range of MAX_MZ_DIFFERENCE), the value will be NaN
# peak_indices[i, j] will contain the index of that peak in 'allPeaks[i]'
peak_masses  = pd.DataFrame(index=reference_df.index, columns=dsSampleCodes,
                            dtype=np.single)
peak_indices = pd.DataFrame(index=reference_df.index, columns=dsSampleCodes,
                            dtype=int)
for s in tqdm(dsSampleCodes, desc='Identifying metabolites'):
    for j, refmass in reference_df['mass'].items():
        diffs = abs(allPeaks[s][:, 0] + proton_mass - refmass)
        peak_idx = np.argmin(diffs)
        if diffs[peak_idx] <= MAX_MZ_DIFFERENCE:
            peak_indices.loc[j, s] = peak_idx
            peak_masses.loc[j, s] = allPeaks[s][peak_idx, 0]
        else:
            peak_indices.loc[j, s] = -1
            peak_masses.loc[j, s] = np.nan

# keep only the reference masses that actually have a 'hit' in at least one
# of the samples, and calculate the median of all samples where a peak was 
# associated with this mass
ref_hits      = (peak_indices != -1).any(1)
peak_indices  = peak_indices.loc[ref_hits, :]
median_masses = peak_masses.loc[ref_hits, :].median(1)
compound_df   = compound_df.loc[ref_hits, :]

#%%
# data_df[i, j] will contain the intensity of the peak which was associated with
# reference mass 'j' in sample 'i'. If there wasn't any close enough mass, 
# we take the median mass of the ions associated with this reference mass
# across all other samples, and find the ion closest to the median (even if 
# it is not actually a peak).

data_df = pd.DataFrame(index=compound_df.index, columns=dsSampleCodes,
                       dtype=np.single)
data_df.index.name = 'ion.code'

for s in tqdm(dsSampleCodes, desc='Creating final matrix'):
    for j, median_mass in median_masses.items():
        peak_idx = peak_indices.loc[j, s]
        if peak_idx != -1:
            # if there is a peak associated with the metabolite,
            # get the intensity of that peak
            data_df.loc[j, s] = allPeaks[s][peak_idx, 1]
        else:
            # otherwise, get the intensity from the closest mz in the raw data
            idx = np.argmin(np.abs(median_mass - dataProfiles[s][:, 0]))
            data_df.loc[j, s] = dataProfiles[s][idx, 1]

#merged = compound_df.join(merged)

#%%
if args.output_path is None:
    args.output_path = os.path.join(os.path.abspath(os.path.curdir),
                                    args.exp_code)

ion_fname = args.output_path + '_ions.csv'
sample_fname = args.output_path + '_samples.csv'
data_fname = args.output_path + '_data.csv'

sys.stderr.write('\nWriting results to output CSV files to path "%s" ... '
                 % args.output_path)
compound_df.to_csv(ion_fname)
data_df.to_csv(data_fname)
sample_df.to_csv(sample_fname)

sys.stderr.write('[DONE]\n')

