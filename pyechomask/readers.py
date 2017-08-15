# -*- coding: utf-8 -*-
"""
.. :module:: readers
    :synopsis: Read in raw files/data objects and output Sv_dict[fq]

            
            Sv_dict[fq]:
                    dict{
                                Sv: numpy array[nSamples by nPings]:          Sv values: (dB re 1m-1)
                       start_depth:            list[floats][nPings]: first sample depth: (m)
                        sample_int:            list[floats][nPings]:    sample interval: (m)
                      pulse_length:            list[floats][nPings]:       pulse length: (ms)
                    }
            
            fq: float: frequency in kHz
            For multi-frequency data, ping number assumed to be equal to column number 
            across all frequencies


| Developed by: Roland Proud (RP) <rp43@st-andrews.ac.uk> 
|               Pelagic Ecology Research Group, University of St Andrews
| Contributors:
|
| Maintained by:
| Modification History:      
|
"""

import gzip
import pickle
import numpy as np


def read_PERGobjs(fp,mask_value = -999):
    '''
    :param fp: filepath to object
    :type  fp: string
    
    parse objects stored in the pelagic ecology research group database
    and output Sv_dict
    
    '''
    
    ## read raw multi-frequency EK60 data
    f   = gzip.open(fp,'rb')
    obj = pickle.load(f,encoding = 'bytes')
    f.close()
    
    Sv_dict = {}
    
    for fq in obj.keys():
        fq = int(fq)
        ## initialize dict
        Sv_dict[fq] = {}
        
        ##get dimensions
        startPing = 0
        maxRow    = 0
        for pingRange in obj[fq].keys():
            endPing         = pingRange[1]
            row_idx,col_idx = obj[fq][pingRange]['idx']
            maxRow          = max(maxRow,np.max(row_idx))
        
        ## initialize grid
        Sv_grid     = np.ones((maxRow + 1,endPing + 1)) * mask_value
        start_depth = []
        samp_int    = []
        tPL         = []
        ## populate grid    
        for pingRange in obj[fq].keys():
            sample_int      = obj[fq][pingRange]['sample_int']
            idx             = obj[fq][pingRange]['idx']
            Sv              = obj[fq][pingRange]['Sv'].astype('f8')
            row_idx,col_idx = idx
            
            ## sv_grid
            Sv_grid[row_idx,col_idx + pingRange[0]] = Sv
            ## depth
            for s in sample_int:
                start_depth = start_depth + list(np.repeat(s[2],(s[1]-s[0]) + 1))
                samp_int    = samp_int + list(np.repeat(s[3],(s[1]-s[0]) + 1))
                tPL         = tPL + list(np.repeat(obj[fq][pingRange]['paras']['pulse_length'],(s[1]-s[0]) + 1))
        
        Sv_dict[fq] = {'Sv':Sv_grid,'start_depth':start_depth,'sample_int':samp_int,'pulse_length':tPL}
        

    return Sv_dict


def remove_background_noise(Sv,mask_value = -999):
    '''
    Add when other readers added
    '''
    
    


