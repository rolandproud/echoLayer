"""
.. :module:: layers
    :synopsis: Contains standardised layer algorithms

              Each method reads raw_reader objects (not yet implemented)
              and dictionary of parameters to produce layer. 
              
              Method format:
 
              def layer-type_unique-name(obj,params):
              '''
              param descs
              ...
              function desc
              '''
              code
              ...
              return layer

             layer-type can be 'binary' (0 or 1), 'flag'(range of ints) 
             or 'cont' (continuous: values range from 0-1)

| Developed by: Roland Proud <rp43@st-andrews.ac.uk> 
|               Pelagic Ecology Research Group, University of St Andrews
| Maintained by:
|
"""

## import packages
import numpy as np

###########################################################################################

def binary_threshold(obj,params):
    '''
    :param obj: echosounder data and parameters
    :type  obj: raw_reader object
    
    :param params: th: threshold-value (dB re 1m^-1)
    :type  params: dictionary of parameter key-value pairs:
                   th: float

    generate threshold mask
    
    '''
    Sv                       = obj[params['fq']]['Sv']  
    layer                    = np.zeros(Sv.shape).astype(int)
    layer[Sv > params['th']] = 1
    return layer


