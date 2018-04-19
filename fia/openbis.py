# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 10:53:11 2015

@author: noore
"""

import psycopg2
import numpy as np
import itertools
from collections import OrderedDict
from tqdm import tqdm
import settings as S

#% -------------------------------------------------------------------------
#% fiaMLconfig - fiaexp database config file
#% -------------------------------------------------------------------------
#% [filedname]_table         database table containing the fields data 
#%                           valid values :  data_set_properties
#%                                           sample_properties
#%
#% [filedname]_pt_id         property type id
#%                           valid values :  all integer
#%
#% [filedname]_customParser  anonymous MATLAB function to parse database
#%                           content
#%                           valid values :  empty
#%                                           valid MATLAB function (@(x)
#%                                           ...)
#%
#% [filedname]_type          data type to cast values into
#%                           valid values :  numeric
#%                                           alphanumeric
#%
#% -------------------------------------------------------------------------
#% Sauerlab setup:
#% -------------------------------------------------------------------------
#% FIELD (openBIS)           TABLE                   PT_ID
#% DS COMMENTS               data_set_properties     2
#% DS SAMPLENAME             sample_properties       1
#% DS SAMPLEPERTURBATION     sample_properties       3
#% DS SAMPLETIME             sample_properties       4
#% DS SAMPLEPOSITION         sample_properties       5
#% DS SAMPLEAMOUNT           sample_properties       6
#% DS SAMPLESPECIES          sample_properties       7
#% DS FILENAME               data_set_properties     5
#% -------------------------------------------------------------------------

# openbis config
STPT_ID = 1
DSTY_ID = 3

def list_to_comma_separated_string(l):
    if type(l[0]) == str:
        return ','.join(map(lambda s: "'" + s + "'", l))
    else:
        return ','.join(map(str, l))

def download_data_profiles(exp_code):
    conn_ob = psycopg2.connect(database=S.OPENBIS_DSN,
                               user=S.OPENBIS_LOGIN,
                               password=S.PASSWORD,
                               host=S.HOST,
                               port=S.PORT)
    conn_mb = psycopg2.connect(database=S.METABOL_DSN,
                               user=S.METABOL_LOGIN,
                               password=S.PASSWORD,
                               host=S.HOST,
                               port=S.PORT)
    cur_ob = conn_ob.cursor()
    cur_mb = conn_mb.cursor()

    # find the experiment ID (an openBIS internal ID, not the code that appears on
    # the website)
    cur_ob.execute("SELECT id FROM experiments WHERE code='%s';" % exp_code)
    exp_id = cur_ob.fetchone()[0]
    
    # get the conversion dictionary from 
    cur_ob.execute("""SELECT   data.code, samples.code 
                      FROM     data, samples
                      WHERE    data.expe_id='%d' AND data.dsty_id=%d
                        AND    data.samp_id = samples.id
                      ORDER BY data.code""" % 
                   (exp_id, DSTY_ID))
    dsCode2smpCode = OrderedDict(cur_ob.fetchall())

    cur_mb.execute("""SELECT   data_sets.perm_id, fia_ms_runs.id
                      FROM     data_sets, fia_ms_runs
                      WHERE    data_sets.perm_id in (%s)
                        AND    data_sets.id = fia_ms_runs.ds_id
                      ORDER BY data_sets.perm_id""" % 
                   list_to_comma_separated_string(list(dsCode2smpCode.keys())))
    dsCode2fiaId = OrderedDict(cur_mb.fetchall())

    if list(dsCode2smpCode.keys()) != list(dsCode2fiaId.keys()):
        raise Exception('Could not find all the datasets in openBIS, aborting...')
    
    dataProfiles = {}
    for dsCode, fia_ms_run_id in tqdm(dsCode2fiaId.items(),
                                      desc='Downloading m/z data'):
        cur_mb.execute("""SELECT   mz, intensities 
                          FROM     fia_profiles
                          WHERE    fia_ms_run_id=%d 
                          ORDER BY low_mz"""
                       % fia_ms_run_id)
        mz, intens = zip(*cur_mb.fetchall())
        mz = map(lambda s: map(float, s.split(',')), mz)
        mz = list(itertools.chain(*mz))
    
        intens = map(lambda s: map(float, s.split(',')), intens)
        intens = list(itertools.chain(*intens))
        profile = np.matrix([mz, intens], dtype=np.single).T
        dataProfiles[dsCode2smpCode[dsCode]] = profile
    
    conn_ob.close()
    conn_mb.close()
    return dataProfiles

def get_sample_names(exp_code):
    conn_ob = psycopg2.connect(database=S.OPENBIS_DSN,
                               user=S.OPENBIS_LOGIN,
                               password=S.PASSWORD,
                               host=S.HOST,
                               port=S.PORT)
    cur_ob = conn_ob.cursor()

    # find the experiment ID (an openBIS internal ID, not the code that appears on
    # the website)
    cur_ob.execute("SELECT id FROM experiments WHERE code='%s';" % exp_code)
    exp_id = cur_ob.fetchone()[0]

    cur_ob.execute("""SELECT   samples.code, sample_properties.value
                      FROM     data, samples, sample_properties
                      WHERE    data.expe_id='%d' AND data.dsty_id=%d
                        AND    data.samp_id = samples.id
                        AND    samples.id = sample_properties.samp_id
                        AND    sample_properties.stpt_id = %d
                      ORDER BY data.code""" % 
                   (exp_id, DSTY_ID, STPT_ID))
    dsCode2name = OrderedDict(cur_ob.fetchall())
    return dsCode2name
    
if __name__ == '__main__':
    dsCode2name = get_sample_names('E222456')
    print(dsCode2name)