#!/usr/bin/python
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
import matplotlib.gridspec as gridspec
sns.set()

parser = argparse.ArgumentParser(description='Analyze Bradford assay (protein quantification) data')
parser.add_argument('excel_file', type=argparse.FileType(mode='r'),
                    help='path to an Excel file containing the Bradford measurements')
parser.add_argument('-o', dest='output_fname', type=str, default=None,
                    help='a filename for writing the output')
parser.add_argument('-s', dest='sunrise', type=bool, default=False,
                    help='use the Sunrise format (i.e. not Tecan Infinite 200M)')
args = parser.parse_args()

# if no output fname is given, use the Excel fname but with a .png extension
if args.output_fname is None:
    output_fname = args.excel_file.name + '.png'
else:
    output_fname = args.output_fname

PROTEIN_CONC_L = 'Protein concentration [$\mu$g/mL]'
STDEV_L = 'Standard dev. [$\mu$g/mL]'
def get_all_tables(fp):
    sheet_to_df = {}
    for sheet_name, sheet_df in pd.read_excel(fp, None, index_col=None, header=None).iteritems():
        first_col = sheet_df[0]
        start_rows = sheet_df[first_col == '<>'].index
        for start_row in start_rows:
            # iterate row letters from the first one after the <> until the
            # end of the table (or stop at 'H').
            end_row = start_row+1
            letter = first_col[end_row]
            while end_row in first_col.index and letter == first_col[end_row]:
                end_row += 1
                letter = chr(ord(letter) + 1)
                
            tbl = sheet_df.iloc[start_row+1:end_row, 1:]
            tbl.columns = sheet_df.iloc[start_row, 1:].fillna(0).apply(int)
            tbl.columns.name = 'col'
            tbl['row'] = first_col[start_row+1:end_row]
            tbl = pd.melt(tbl, id_vars=['row'])
            tbl['well'] = tbl['row'] + map(str, tbl['col'])
            tbl = tbl[~pd.isnull(tbl['value'])]
            tbl.index = tbl['well']
            tbl = tbl[['value']]
            sheet_to_df.setdefault(sheet_name, []).append(tbl)
    return sheet_to_df

sheet_to_df = get_all_tables(args.excel_file)
if args.sunrise:
    df_595, df_label, df_dilution = sheet_to_df.values()[0]
    df_ratio = df_595
    ylabel = 'Absorbance (595nm)'
else:
    if 'Sheet2' not in sheet_to_df or len(sheet_to_df['Sheet2']) != 2:
        raise ValueError('The Tecan measurements must be found in "Sheet2", for both 590nm and 450nm')
    if 'Sheet1' not in sheet_to_df:
        raise ValueError('The labels and dilution ratios must be found in "Sheet1"')
    
    df_590, df_450 = sheet_to_df['Sheet2']
    df_label, df_dilution = sheet_to_df['Sheet1']

    df_ratio = df_590/df_450
    ylabel = 'Absorbance ratio (590nm / 450nm)'

df_ratio.rename(columns={'value': 'ratio'}, inplace=True)
df_label.rename(columns={'value': 'label'}, inplace=True)
df_dilution.rename(columns={'value': 'dilution'}, inplace=True)

df_data = df_label.join(df_dilution).join(df_ratio).sort_index()

# find all the "standard" samples, and create the calibration curve
df_std = df_data[df_data['label'] == 'STD']
if df_std.shape[0] == 0:
    raise ValueError('All standards must be marked by the same label: "STD"')
try:
    slope, intercept, _, _, _ = linregress(df_std['dilution'].tolist(), df_std['ratio'].tolist())
except AttributeError as e:
    print df_std
    raise e
conc_to_ratio = lambda c : c*slope + intercept
ratio_to_conc = lambda r : (r - intercept) / slope
df_data['conc'] = df_data['ratio'].apply(ratio_to_conc)

#%%
fig = plt.figure(figsize=(10, 4)) 
gs = gridspec.GridSpec(2, 3, width_ratios=[3.5, 0.3, 6], height_ratios=[3, 0.1], top=0.9, bottom=0.2)
#fig, axs = plt.subplots(1, 2, figsize=(15, 5))
ax0 = plt.subplot(gs[0,0])
ax0.plot(df_std['dilution'] + 1e-2, df_std['ratio'], 'o')
ax0.plot(np.linspace(1e-2, 1600, 100), map(conc_to_ratio, np.linspace(1e-2, 1600, 100)), '-')
ax0.set_xlabel(PROTEIN_CONC_L)
ax0.set_ylabel(ylabel)
ax0.set_xlim(0, None)
ax0.set_ylim(0, None)

df_samples = df_data[df_data['label'] != 'STD']
df_samples.loc[:, 'conc'] = df_samples['conc'] * df_samples['dilution']

ax1 = plt.subplot(gs[:,2])
res_df = pd.DataFrame(index=df_samples['label'].unique(), columns=[PROTEIN_CONC_L, STDEV_L])
res_df.index.name = 'Sample'
conc_groups = df_samples.groupby('label')['conc']
res_df[PROTEIN_CONC_L] = conc_groups.sum() / conc_groups.count()
for l in res_df.index:
    conc_vec = df_samples.loc[df_samples['label'] == l, 'conc']
    res_df.loc[l, STDEV_L] = np.round(conc_vec.std(), 1)
res_df[PROTEIN_CONC_L] = res_df[PROTEIN_CONC_L].round(1)
res_df[PROTEIN_CONC_L].plot(kind='bar', yerr=res_df[STDEV_L], 
                            ax=ax1, color=(0.9, 0.6, 0.65), table=res_df.T,
                            linewidth=0)
ax1.set_xlabel('')
ax1.set_ylabel(PROTEIN_CONC_L)
ax1.get_xaxis().set_visible(False)
fig.savefig(output_fname, pgi=300)
