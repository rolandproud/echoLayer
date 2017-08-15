# -*- coding: utf-8 -*-
"""
.. :module:: masks
    :synopsis: Contains standardised mask methods

              Each method reads Sv_dict objects (see readers.py)
              and dictionary of parameters to produce a mask. 
              
              Method format:
 
                def [mask-type]_[unique-name](Sv_dict,params):
                    '''
                    :param Sv_dict: contains Sv grid, depth
                    :type  Sv_dict: dictionary
                    
                    :param params: fq: frequency (kHz)
                    :type  params: dictionary of parameter key-value pairs:
                                   fq: float
                                   ...
                
                    [description]
                    
                    defined by [initials of developer]
                    
                    status: [status(dev,test,product)]
                    
                    '''
                    ## get Sv grid and create mask grid
                    Sv   = Sv_dict[params['fq']]['Sv'] 
                    mask = np.zeros(Sv.shape).astype(int)
                    
                    ## mask code   
                    
                    
                    code...
                    
                    
                    ## add mask to Sv_dict
                    Sv_dict[params['fq']]['mask'] = mask 
                    
                    return Sv_dict[params['fq']]
                
             

             mask-type can be 'binary' (0 or 1), 'flag'(range of ints) 
             or 'cont' (continuous: values range from 0-1)
             
             for binary masks, 1 = signal; 0 = noise

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

################################################################## background noise

## Add function to reader module
## background noise is already removed from PERGobjects

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
    
    status: product
    
    '''
    ## get Sv grid and create mask grid
    Sv                            = Sv_dict[params['fq']]['Sv'] 
    mask                          = np.zeros(Sv.shape).astype(int)
    
    ## mask code
    mask[Sv > params['th']]       = 1
    
    ## add mask to Sv_dict
    Sv_dict[params['fq']]['mask'] = mask 
        
    return Sv_dict[params['fq']]

## detect aggregates/SSLs

################################################################### noise masks

## transmit pulse and near-field

def binary_pulse(Sv_dict,params):
    '''
    :param Sv_dict: contains Sv grid, depth
    :type  Sv_dict: dictionary
    
    :param params: fq: frequency (kHz)
    :type  params: dictionary of parameter key-value pairs:
                   fq: float
                   ...

    generate pulse mask, mask pulse and surface noise
    
    defined by RP
    
    status: dev
    
    '''
    ## get Sv grid and create mask grid
    Sv   = Sv_dict[params['fq']]['Sv'] 
    mask = np.ones(Sv.shape).astype(int)
    
    ## mask code   
    samples,pings = Sv.shape
    for p in range(pings):
        idx           = np.where(Sv[:,p] == params['mask_value'])[0][0]
        mask[0:idx,p] = 0
        
    ## add mask to Sv_dict
    Sv_dict[params['fq']]['mask'] = mask 
    
    return Sv_dict[params['fq']]


## surface (bubbles/airation)

## seabed

## false-bottom

## dropped ping/attenuated signal

## transient noise/noise spike

## impulse/interference - regular discrete pulses of sound from external source

## lowered instrument (CTD etc.)







































