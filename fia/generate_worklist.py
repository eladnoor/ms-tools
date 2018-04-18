#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 10:41:42 2018

@author: noore
"""
import pandas as pd
import argparse
import os

DATA_PATH = 'D:/MassHunter/data/'
METHOD_PATH = 'D:/MassHunter/data/all/method repository/fia_ctc_neg_AF_double.m'
SAMPLE_RUNTIME_MIN = 1.56


def GetArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--trays', action='append',
                        help='list of tray numbers')
    parser.add_argument('-w', '--water', type=bool,
                        help='water plate in tray 1')
    parser.add_argument('--sbSkipCheck', action='store_true',
                        help='skip check tray numbers')
    parser.add_argument('-u', '--username', type=str, required=True,
                        help='username is only used for determining output directory')
    parser.add_argument('-d', '--datasetname', type=str, required=True,
                        help='datasetname is only used for determining output directory')
    parser.add_argument('input_fname', type=argparse.FileType('r'),
                        help='CSV input file')
    args = parser.parse_args()
    return args

def ParseTrays(args):
    if args.trays:
        trays = set(map(int, args.trays))
    else:
        trays = set()
    if args.water:
        trays.update(1)
    return trays

def ParseSamplesFile(args):
    """
        Read the sample table from the CSV input file
    """
    def FindColumn(df, col_names, raise_exception=True):
        for name in col_names:
            if name in df.columns:
                return name
        if raise_exception:
            raise ValueError('A column named "%s" must be present in the Excel file'
                             % col_names[0])
        else:
            return None

    input_df = pd.read_csv(args.input_fname)
    
    # obligatory columns
    posCol = FindColumn(input_df,
                        ['position', 'pos', 'well'])
    nameCol = FindColumn(input_df,
                         ['name', 'samplename', 'sample', 'sample name'])
    pertCol = FindColumn(input_df,
                         ['perturbation', 'condition', 'mut'])
    input_df.rename(columns={posCol: 'position', nameCol: 'name', pertCol:'perturbation'},
                     inplace=True)

    if not input_df.name.is_unique:
        raise ValueError('The "name" column must only contain unique values')

    return input_df

#    # optional columns
#    plateCol = FindColumn(raw_input, ['plate'], raise_exception=False)
#    groupCol = FindColumn(raw_input, ['group'], raise_exception=False)
#    amountCol = FindColumn(raw_input,
#                           ['amount', 'od', 'weight'], raise_exception=False)
#    grCol = FindColumn(raw_input,
#                       ['gr', 'growth', 'rate', 'time'], raise_exception=False)
#    speciesCol = FindColumn(raw_input,
#                            ['species', 'sample species'], raise_exception=False)


###############################################################################
if __name__ == '__main__':
    args = GetArgs()
    trays = ParseTrays(args)
    input_df = ParseSamplesFile(args)
    
    if 'plate' in input_df:
        n_plates = input_df.plate.nunique()
        if args.sbSkipCheck:
            trays = set(range(1, n_plates+1))
        elif len(trays) < n_plates:
            raise ValueError('You need to select at least %d trays!' % n_plates)
    else:
        input_df['plate'] = 1
    
    path = os.path.join(DATA_PATH, args.username, args.datasetname)
    prefix = args.datasetname

    # create a samples DataFrame (for OpenBIS)
    samples_df = pd.DataFrame(index=input_df.index)
    samples_df['SAMPLEPOSITION'] = input_df.position.str.cat(input_df.plate.map(str), delimiter='-')
