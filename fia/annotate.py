#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 23:54:03 2015

@author: eladn
"""

import numpy as np
import sys
import csv
import os
import argparse
from progressbar import ProgressBar
from openbis import download_data_profiles

# script parameters (as determined by Mattia)
MIN_PEAK_SIZE = 5000
MAX_MZ_DIFFERENCE = 0.003
MAX_REF_MASS = 1000
base_path = os.path.split(os.path.realpath(__file__))[0]
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

parser = argparse.ArgumentParser(description='Download FIA raw data from openBIS')
parser.add_argument('exp_code', type=str,
                    help='the openBIS experiment ID')
parser.add_argument('-o', dest='output_fname', type=str, default=None,
                    help='a filename for writing the output')
args = parser.parse_args()

dsCode, dsSampleCode, dataProfiles = download_data_profiles(args.exp_code)
n_samples = len(dataProfiles)

# all_peaks is a list of matrices (one per sample) of the ion mz and intensity
# only for the peaks (i.e. local maxima)
all_peaks = []

#%% identify peaks (local maxima)
sys.stderr.write('\nCentroids identification\n')
with ProgressBar(max_value=n_samples) as progress:
    for i, dataset in enumerate(dataProfiles):
        # find all the values that are local maxima and pass the threshold
        idxs = findpeaks(dataset[:, 1])
        idxs = filter(lambda j : dataset[j, 1] >= MIN_PEAK_SIZE, idxs)
        all_peaks.append(dataset[idxs, :])
        progress.update(i)

#%% Use the reference table to associate peaks to compounds, by minimum mass distance
sys.stderr.write('\nMetabolites identification\n')

with open(REFERENCE_MASS_FNAME, 'r') as fp:
    # we didn't measure masses above 1000 Da
    reference_mass_dlist = filter(lambda x:float(x['mass']) < MAX_REF_MASS, csv.DictReader(fp))
    reference_masses = np.array(map(lambda x: float(x['mass']), reference_mass_dlist))
    reference_masses -= reference_masses[0] # subtract the mass of H+ (i.e. look for deprotonated masses)

# peak_masses[i, j] will contain the exact mass of the peak which is closest
# to reference 'j' in sample 'i'. If there is no peak which is close enough
# (i.e. in the range of MAX_MZ_DIFFERENCE), the value will be NaN
# peak_indices[i, j] will contain the index of that peak in 'all_peaks[i]'
peak_masses = np.zeros((n_samples, len(reference_mass_dlist)), dtype=np.single)
peak_indices = np.zeros((n_samples, len(reference_mass_dlist)), dtype=int)
with ProgressBar(max_value=n_samples) as progress:
    for i in xrange(n_samples):
        for j, refmass in enumerate(reference_masses):
            diffs = abs(all_peaks[i][:, 0] - refmass)
            peak_idx = np.argmin(diffs)
            if diffs[peak_idx] <= MAX_MZ_DIFFERENCE:
                peak_indices[i, j] = peak_idx
                peak_masses[i, j] = all_peaks[i][peak_idx, 0]
            else:
                peak_indices[i, j] = -1
                peak_masses[i, j] = np.nan
        progress.update(i)

# keep only the reference masses that actually have a 'hit' in at least one
# of the samples.
ref_hits             = list(np.where(np.any(peak_indices != -1, 0))[0])
median_masses        = np.nanmedian(peak_masses[:, ref_hits], axis=0)
peak_indices         = peak_indices[:, ref_hits]
reference_mass_dlist = [reference_mass_dlist[i] for i in ref_hits]
n_references         = len(ref_hits)

#%%
sys.stderr.write('\nCreating final matrix\n')

# merged[i, j] will contain the intensity of the peak which was associated with
# reference mass 'j' in sample 'i'. If there wasn't any close enough mass, 
# we take the median mass of the ions associated with this reference mass
# across all other samples, and find the ion closest to the median (even if 
# it is not actually a peak).

merged = np.zeros((n_samples, n_references), dtype=np.single)
with ProgressBar(max_value=n_samples) as progress:
    for i, dataset in enumerate(dataProfiles): # 13800 sec
        for j, median_mass in enumerate(median_masses):
            peak_idx = peak_indices[i, j]
            if peak_idx != -1:
                # if there is a peak associated with the metabolite,
                # get the intensity of that peak
                merged[i, j] = all_peaks[i][peak_idx, 1]
            else:
                # otherwise, get the intensity from the closest mz in the raw data
                idx = np.argmin(np.abs(median_mass - dataset[:, 0]))
                merged[i, j] = dataset[idx, 1]
        progress.update(i)

#%%
if args.output_fname is None:
    args.output_fname = args.exp_code + '.csv'

sys.stderr.write('\nWriting results to output CSV file "%s" ... ' % args.output_fname)
with open(args.output_fname, 'w') as fp:
    csv_output = csv.writer(fp)
    ref_keys = ['name', 'keggId', 'formula', 'mass']
    csv_output.writerow(ref_keys + map(str, xrange(1, n_samples+1)))
    for j, rdict in enumerate(reference_mass_dlist):
        csv_output.writerow(map(rdict.get, ref_keys) + list(merged[:, j].flat))

sys.stderr.write('[DONE]\n')
