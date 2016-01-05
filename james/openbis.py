# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 10:53:11 2015

@author: noore
"""

import psycopg2
import pandas as pd

# openbis config
HOST          = 'baobab.ethz.ch'        # hostname of the server running openBIS
PORT          = 5432                    # port for openBIS database
PASSWORD      = 'openbisyeastx'         # password for openBIS
OPENBIS_DSN   = 'openbis_productive'    # data source name of the openBIS database
OPENBIS_LOGIN = 'openbis_readonly'      # username for openBIS
DSTY_ID       = 4
STPT_DICT     = {1: 'Name',
                 2: 'Species/Strain/Cell',
                 3: 'Perturbation',
                 4: 'Time point (min)',
                 5: 'Position',
                 6: 'Amount'}

PERTURBATION_COL_NAME = STPT_DICT[3]
TIME_COL_NAME = STPT_DICT[4]
STRAIN_COL_NAME = STPT_DICT[2]

FLOAT_COLUMNS = [4, 6]

def download_metadata(exp_code):
    conn_ob = psycopg2.connect(database=OPENBIS_DSN,
                               user=OPENBIS_LOGIN,
                               password=PASSWORD,
                               host=HOST,
                               port=PORT)
    
    sql = """SELECT   d.samp_id samp_id, d.code perm_id, d.id ds_id, s.code, sp.stpt_id, sp.value
             FROM     experiments e, data d, samples s, sample_properties sp
             WHERE    e.code='%s'
                AND   d.expe_id = e.id
                AND   d.dsty_id = %d
                AND   s.id = d.samp_id
                AND   s.id = sp.samp_id
          """ % (exp_code, DSTY_ID)
    result_df = pd.read_sql_query(sql, conn_ob)
    conn_ob.close()
    
    # reshape the result so that the perm_id is the index, and the columns are
    # the different property types
    new_result_df = pd.DataFrame(index=sorted(result_df['perm_id'].unique()))
    
    for i in sorted(result_df['stpt_id'].unique()):
        tmp = result_df[result_df['stpt_id'] == i].drop_duplicates()
        new_result_df[STPT_DICT[i]] = tmp.set_index('perm_id')['value']
    new_result_df.sort_index(inplace=True)
    
    for i in FLOAT_COLUMNS:
        if STPT_DICT[i] in new_result_df.columns:
            new_result_df.loc[:, STPT_DICT[i]] = \
                new_result_df.loc[:, STPT_DICT[i]].apply(float)

    return new_result_df
    
if __name__ == '__main__':
    # example
    samples = download_metadata('E221350')
    print samples