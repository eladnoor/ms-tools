# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 10:53:11 2015

@author: noore
"""

import psycopg2
import pandas as pd

# openbis config
host          = 'baobab.ethz.ch'        # hostname of the server running openBIS
port          = 5432                    # port for openBIS database
driver        = 'org.postgresql.Driver' # driver used to access openBIS
password      = 'openbisyeastx'         # password for openBIS
openbis_dsn   = 'openbis_productive'    # data source name of the openBIS database
openbis_login = 'openbis_readonly'      # username for openBIS
dsty_id       = 4

stpt_dict     = {1: 'Name',
                 2: 'Species/Strain/Cell',
                 3: 'Perturbation',
                 4: 'Time point (min)',
                 5: 'Position',
                 6: 'Amount'}

def download_metadata(exp_code):
    conn_ob = psycopg2.connect(database=openbis_dsn,
                               user=openbis_login,
                               password=password,
                               host=host,
                               port=port)
    
    sql = """SELECT   d.samp_id samp_id, d.code perm_id, d.id ds_id, s.code, sp.stpt_id, sp.value
             FROM     experiments e, data d, samples s, sample_properties sp
             WHERE    e.code='%s'
                AND   d.expe_id = e.id
                AND   d.dsty_id = 4
                AND   s.id = d.samp_id
                AND   s.id = sp.samp_id
          """ % (exp_code)
    result_df = pd.read_sql_query(sql, conn_ob)
    conn_ob.close()
    
    # reshape the result so that the perm_id is the index, and the columns are
    # the different property types
    new_result_df = pd.DataFrame(index=sorted(result_df['perm_id'].unique()))
    
    for i in sorted(result_df['stpt_id'].unique()):
        tmp = result_df[result_df['stpt_id'] == i].drop_duplicates()
        new_result_df[stpt_dict[i]] = tmp.set_index('perm_id')['value']
    new_result_df.sort_index(inplace=True)
    
    return new_result_df
    
if __name__ == '__main__':
    samples = download_metadata('E221238')