"""
.. :module:: layers
    :synopsis: Contains standardised layer algorithms

              Each method reads raw_reader objects (not yet implemented)
              and dictionary of parameters to produce a layer. 
              
              Method format:
 
              def layer-type_unique-name(obj,params):
              '''
              param descs
              ...
              function desc
              '''
              code
              ...
              return {'layer':layer,'ping_time':ping_time,\
                'ping_sample_time':ping_sample_time,'pulse_length':pulse_length}
                
             layer: 2D numpy.array contains layer values
             from raw_reader object:
             ping_time: time of each ping,list of datetime.datetime objects 
             ping_sample_time: 2D numpy.array, sample time relative to ping  
             pulse_lenth: pulse length from calibration parameters

             layer-type can be 'binary' (0 or 1), 'flag'(range of ints) 
             or 'cont' (continuous: values range from 0-1)

| Developed by: Roland Proud <rp43@st-andrews.ac.uk> 
|               Pelagic Ecology Research Group, University of St Andrews
| Maintained by:
| Modification History:
| 07/07/2017   Updated standard method format to include ping/sample time and 
|              pulse length (needed to merge layers of different grid sizes).       
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
    ##Define Layer
    Sv                       = obj[params['fq']]['Sv']  
    layer                    = np.zeros(Sv.shape).astype(int)
    layer[Sv > params['th']] = 1
    ## Define ping time, ping sample time and pulse length (from raw_reader obj)
    ping_time        = obj[params['fq']]['dts']
    ping_sample_time = obj[params['fq']]['sample_time']
    pulse_length     = obj[params['fq']]['pulse_length']
    
    return {'layer':layer,'ping_time':ping_time,\
    'ping_sample_time':ping_sample_time,'pulse_length':pulse_length}


