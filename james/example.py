# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 18:45:38 2015

@author: noore
"""
import argparse
from openbis import download_metadata

parser = argparse.ArgumentParser(description='Download metadata of all samples associated with experiment')
parser.add_argument('exp_code', type=str,
                    help='the openBIS experiment ID')
parser.add_argument('-o', dest='output_fname', type=str, default=None,
                    help='a filename for writing the output')
args = parser.parse_args()

samples_df = download_metadata(args.exp_code)
samples_df.rename(columns={'perm_id': 'sample', 'Perturbation': 'exp', 'Time point (min)': 'time (min)'}, inplace=True)
samples_df.to_csv(args.output_fname)