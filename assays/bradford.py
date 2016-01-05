# -*- coding: utf-8 -*-
"""
Created on Tue Jan  5 10:42:05 2016

@author: noore
"""

import argparse
import pandas as pd
import numpy as np
from scipy.stats import linregress
import seaborn as sns
import matplotlib.pyplot as plt
sns.set()

parser = argparse.ArgumentParser(description='Analyze Bradford assay (protein quantification) data')
parser.add_argument('excel_file', type=argparse.FileType(mode='r'),
                    help='path to an Excel file containing the Bradford measurements')
parser.add_argument('output_fname', type=argparse.FileType(mode='w'),
                    help='path for writing the figure (SVG format)')
args = parser.parse_args()

def get_all_tables(fp):
    sheet_to_df = {}
    for sheet_name, sheet_df in pd.read_excel(fp, None, index_col=None, header=None).iteritems():
        start_rows = sheet_df[sheet_df[0] == '<>'].index
        for start_row in start_rows:
            # iterate all the letters from A to H and stop when the value doesn't match
            # the expected row name
            for i, end_row in enumerate(range(start_row+1, sheet_df.shape[0])):
                if (sheet_df.iloc[end_row, 0] != chr(ord('A') + i)):
                    break
            table = sheet_df.iloc[start_row+1:end_row, 1:]
            table.columns = sheet_df.iloc[start_row, 1:].apply(int)
            table.columns.name = 'col'
            table['row'] = sheet_df.iloc[start_row+1:end_row, 0]
            table = pd.melt(table, id_vars=['row'])
            table['well'] = table['row'] + map(str, table['col'])
            table = table[~pd.isnull(table['value'])]
            table.index = table['well']
            table = table[['value']]
            sheet_to_df.setdefault(sheet_name, []).append(table)
    return sheet_to_df

sheet_to_df = get_all_tables(args.excel_file)
if 'Sheet2' not in sheet_to_df or len(sheet_to_df['Sheet2']) != 2:
    raise ValueError('The Tecan measurements must be found in "Sheet2", for both 590nm and 450nm')
if 'Sheet1' not in sheet_to_df:
    raise ValueError('The labels and dilution ratios must be found in "Sheet1"')

df_590, df_450 = sheet_to_df['Sheet2']
df_label, df_dilution = sheet_to_df['Sheet1']

df_ratio = df_590/df_450
df_ratio.rename(columns={'value': 'ratio'}, inplace=True)
df_label.rename(columns={'value': 'label'}, inplace=True)
df_dilution.rename(columns={'value': 'dilution'}, inplace=True)

df_data = df_label.join(df_dilution).join(df_ratio).sort_index()

# find all the "standard" samples, and create the calibration curve
df_std = df_data[df_data['label'] == 'STD']
slope, intercept, _, _, _ = linregress(df_std['dilution'], df_std['ratio'])
conc_to_ratio = lambda c : c*slope + intercept
ratio_to_conc = lambda r : (r - intercept) / slope
df_data['conc'] = df_data['ratio'].apply(ratio_to_conc)

fig, axs = plt.subplots(1, 2, figsize=(8, 4))
axs[0].plot(df_std['dilution'] + 1e-2, df_std['ratio'], 'o')
axs[0].plot(np.linspace(1e-2, 1600, 100), map(conc_to_ratio, np.linspace(1e-2, 1600, 100)), '-')
axs[0].set_xlabel('Protein concentration [$\mu$g/mL]')
axs[0].set_ylabel('Absorbance ratio (590nm / 450nm)')

df_samples = df_data[df_data['label'] != 'STD']
df_samples.loc[:, 'conc'] = df_samples['conc'] * df_samples['dilution']
df_samples.pivot(columns='label', values='conc').plot(kind='box', ax=axs[1], rot=45)
axs[1].set_ylabel('Protein concentration [$\mu$g/mL]')
fig.tight_layout()
fig.savefig(args.output_fname, format='svg')