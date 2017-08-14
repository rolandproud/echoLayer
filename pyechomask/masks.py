# -*- coding: utf-8 -*-
"""
.. :module:: masks
    :synopsis: Contains standardised mask methods

              Each method reads Sv_dict objects (see readers.py)
              and dictionary of parameters to produce a mask. 
              
              Method format:
 
              def mask-type_unique-name(Sv_dict,params):
              '''
              param descs
              ...
              function desc
              ...
              defined by: initials of developer
              '''
              code
              ...
              return Sv_dict[fq]
                
             

             mask-type can be 'binary' (0 or 1), 'flag'(range of ints) 
             or 'cont' (continuous: values range from 0-1)

| Developed by: Roland Proud (RP) <rp43@st-andrews.ac.uk> 
|               Pelagic Ecology Research Group, University of St Andrews
| Contributors:
|
| Maintained by:
| Modification History:      
|
"""

## import packages
import numpy as np

################################################################## signal masks

def binary_threshold(Sv_dict,params):
    '''
    :param Sv_dict: contains Sv grid, depth
    :type  Sv_dict: dictionary
    
    :param params: th: threshold-value (dB re 1m^-1)
    :type  params: dictionary of parameter key-value pairs:
                   th: float

    generate threshold mask
    
    defined by RP
    
    '''
    ##Define Layer
    Sv                            = Sv_dict[params['fq']]['Sv'] 
    mask                          = np.zeros(Sv.shape).astype(int)
    mask[Sv > params['th']]       = 1
    Sv_dict[params['fq']]['mask'] = mask 
        
    return Sv_dict[params['fq']]

## detect aggregates/SSLs

################################################################### noise masks
####NOTE Background noise already removed.



## transmit pulse and near-field

## surface noise (bubbles/airation)

## seabed

## false-bottom

## dropped ping/attenuated signal

## transient noise/noise spike

## impulse/interference - regular discrete pulses of sound from external source

## lowered instrument (CTD etc.)







































