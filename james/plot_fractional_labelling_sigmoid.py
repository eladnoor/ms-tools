# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 14:25:17 2015

@author: noore
"""
import pandas as pd
import numpy as np
from isotope_util import compute_fractions
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

sns.set()
sns.set_style("whitegrid")

# CURVE_P0 = (1.0, 1.0, 1.0, 1.0)
# def sigmoid(x, a0, a1, a2, a3):
#  return (a0 + a1 * x) / (1.0 + a2 * np.exp(-a3 * x))

CURVE_A0 = np.array([1.0, 1.0, 1.0, 1.0, 1.0], dtype=float)


def sigmoid(x, a0, a1, a2, a3, a4):
    return a0 * (1.0 - np.exp(a1 - a2 * x)) / (1.0 + np.exp(a3 - a4 * x))

compounds, maxisos = zip(*[('G6P', 6), ('F6P', 6), ('FBP', 6), ('DHAP', 3), ('xPG', 3), ('PEP', 3)])
compound2maxiso = dict(zip(compounds, maxisos))

colors = ["dusty purple", "windows blue", "teal green", "scarlet", "purplish pink", "orange"]
colpalette = sns.xkcd_palette(colors)
# colpalette = sns.color_palette("PuBuGn_d", len(compounds))

exps = ['eca', 'bsa', 'wt']
titles = [r'$\Delta$pfkA$\Delta$pfkB  +  pfkA from E. coli',
          r'$\Delta$pfkA$\Delta$pfkB  +  pfkA from B. subtilis',
          'Wild-type E. coli']
          
count_df = pd.DataFrame.from_csv('integration_results.txt', sep='\t', index_col=0)
count_df.fillna(0, inplace=True)

samples_df = pd.DataFrame.from_csv('samples.csv', index_col=0)

# remove time point - 45 min, seems to be a mistake
count_df = count_df.loc[samples_df['time (min)'] != 45, :]

for c, maxiso in compound2maxiso.iteritems():
    for i in xrange(maxiso+1):
        col_name = '%s%d' % (c, i)
        if col_name not in count_df.columns:
            count_df[col_name] = 0

# correct the relative abundances for the natural abundance of 13C atoms
corrected_df = count_df * 0
for c, maxiso in compound2maxiso.iteritems():
    cols = ['%s%d' % (c, i) for i in xrange(maxiso+1)]
    for ind in count_df.index:
        counts = [count_df.loc[ind, col] for col in cols]
        fractions = np.array(compute_fractions(counts))
        fractions[fractions < 0] = 0
        fractions /= fractions.sum()
        corrected_df.loc[ind, cols] = fractions

# %% calculate the fractional labeling of each compound (total 13C out of all C)
fraction_df = pd.DataFrame(index=count_df.index, columns=compounds, dtype=float)
for c in compounds:
    cols = ['%s%d' % (c, i) for i in xrange(compound2maxiso[c] + 1)]
    mat = np.matrix(xrange(compound2maxiso[c] + 1)).T
    fraction_df.loc[:, c] = corrected_df.loc[:, cols].dot(mat) / compound2maxiso[c]

fraction_df = fraction_df.join(samples_df)

# %% calculate means and standard deviations
times = list(sorted(set(fraction_df['time (min)'])))

# %% plot time course of all metabolites
fig = plt.figure(figsize=(15, 4))
for i, exp in enumerate(exps):
    ax = fig.add_subplot(1, len(exps), i + 1)
    exp_df = fraction_df[fraction_df['exp'] == exp]

    # calculate mean and std for each time series
    fraction_mean_df = pd.DataFrame(index=times, columns=compounds, dtype=float)
    fraction_std_df = pd.DataFrame(index=times, columns=compounds, dtype=float)
    for t in times:
        idx = (exp_df['time (min)'] == t)
        fraction_mean_df.loc[t, :] = exp_df.loc[idx, :].mean(0)
        fraction_std_df.loc[t, :] = exp_df.loc[idx, :].std(0)

    exp_df.plot(x='time (min)', style='.', ax=ax, color=colpalette, legend=False)

    # curve-fitting for measured data using sigmoid functions
    all_times = np.linspace(min(times), max(times)*1.05, 300)
    compound2popt = {}
    for j, cpd in enumerate(compounds):
        xvalues = times
        yvalues = list(fraction_mean_df.loc[times, cpd])
        try:
            popt, pcov = curve_fit(sigmoid, xvalues, yvalues, p0=CURVE_A0, maxfev=10000)
            compound2popt[cpd] = popt
            f = lambda t: sigmoid(t, *popt)
            if cpd == 'PEP':
                linewidth = 2
                alpha = 1
            else:
                linewidth = 1.5
                alpha = 0.5
            ax.plot(all_times, map(f, all_times), '-', linewidth=linewidth, alpha=alpha,
                    color=colpalette[j])
        except RuntimeError:
            compound2popt[cpd] = CURVE_A0

    final_values = [sigmoid(max(times), *compound2popt[cpd]) for cpd in compounds]
    for j, cpd in enumerate(compounds):
        y = final_values[j]
        ax.text(max(times)*1.1, y-0.03, cpd, color=colpalette[j], fontsize=10)

    ax.set_xlabel('time (min)')
    ax.set_ylabel('fractional labelling')
    ax.xaxis.grid(False)
    ax.set_ylim(0, 0.7)
    ax.set_xlim(min(all_times)-10, max(all_times) + 20)
    ax.set_title(titles[i])
fig.tight_layout()
fig.savefig('labeling.svg')
fig.savefig('labeling.pdf')
