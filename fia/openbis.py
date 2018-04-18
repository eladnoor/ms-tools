# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 10:53:11 2015

@author: noore
"""

import psycopg2
import types
import sys
import numpy as np
import itertools
from collections import OrderedDict
from tqdm import tqdm

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
host      = 'baobab.ethz.ch'  # hostname of the server running openBIS
port      = 5432
driver    = 'org.postgresql.Driver'                  # driver used to access openBIS
use_ssl   = False                                    # Boolean specifying the use of SSl
password  = 'openbisyeastx'                          # password for openBIS
openbis_dsn   = 'openbis_productive'                     # data source name of the openBIS database
openbis_login = 'openbis_readonly'                       # username for openBIS
metabol_dsn   = 'metabol_productive'                     # data source name of the metabol database
metabol_login = 'metabol_readonly'                       # username for openBIS
exty_id = 1
saty_id = 1
dsty_id = 3

#% DS COMMENTs / P[xxx]#[xxx][a-z]
dsPlateStr_table = 'data_set_properties'
dsPlateStr_pt_id = 2

#% DS PERTURBATION
dsSamplePerturbation_table = 'sample_properties'
dsSamplePerturbation_pt_id = 3
dsSamplePerturbation_customParser = ''
dsSamplePerturbation_type = 'alphanumeric'

#% DS SAMPLENAME
dsSampleName_table = 'sample_properties'
dsSampleName_pt_id = 1
dsSampleName_customParser = ''
dsSampleName_type = 'alphanumeric'

#% DS SAMPLEAMOUNT
dsSampleAmount_table = 'sample_properties'
dsSampleAmount_pt_id = 6
dsSampleAmount_customParser = ''
dsSampleAmount_type = 'numeric'

#% DS CONTROL
dsControl_table = 'sample_properties'
dsControl_pt_id = 7
dsControl_customParser = '@parseControls'
dsControl_type = 'numeric'

#% DS USER1
dsUser1_table = 'sample_properties'
dsUser1_pt_id = 1
dsUser1_customParser = ''
dsUser1_type = 'alphanumeric'

#% DS USER2
dsUser2_table = 'sample_properties'
dsUser2_pt_id = 4
dsUser2_customParser = ''
dsUser2_type = 'numeric'

#% DS USER3
dsUser3_table = 'sample_properties'
dsUser3_pt_id = 6
dsUser3_customParser = ''
dsUser3_type = 'numeric'

#% DS USER4
dsUser4_table = 'data_set_properties'
dsUser4_pt_id = 2
dsUser4_customParser = ''
dsUser4_type = 'alphanumeric'

def list_to_comma_separated_string(l):
    if type(l[0]) == str:
        return ','.join(map(lambda s: "'" + s + "'", l))
    else:
        return ','.join(map(str, l))

def download_data_profiles(exp_code):
    conn_ob = psycopg2.connect(database=openbis_dsn,
                               user=openbis_login,
                               password=password,
                               host=host,
                               port=port)
    conn_mb = psycopg2.connect(database=metabol_dsn,
                               user=metabol_login,
                               password=password,
                               host=host,
                               port=port)
    cur_ob = conn_ob.cursor()
    cur_mb = conn_mb.cursor()

    # find the experiment ID (an openBIS internal ID, not the code that appears on
    # the website)
    cur_ob.execute("SELECT id FROM experiments WHERE code='%s';" % exp_code)
    exp_id = cur_ob.fetchone()[0]
    
    sys.stderr.write('\nDownloading datasets associated to experiment %s\n' % exp_code)
    
    # get the conversion dictionary from 
    cur_ob.execute("""SELECT   data.code, samples.code 
                      FROM     data, samples
                      WHERE    data.expe_id='%d' AND data.dsty_id=%d
                        AND    data.samp_id = samples.id
                      ORDER BY data.code""" % 
                   (exp_id, dsty_id))
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
    for i, (dsCode, fia_ms_run_id) in enumerate(tqdm(dsCode2fiaId.items())):
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