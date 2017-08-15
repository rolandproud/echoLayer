'''

Examples of using Mask class and methods

Modification History:

'''

## import packages
import os
import pickle    
import matplotlib.pyplot as plt
import glob
import gzip
import numpy as np

## import user modules
import pyechomask
from pyechomask.mask import Mask
from pyechomask.masks import *
from pyechomask.readers import *
from pyechomask.plotting import *

## create instance of mask class
mask_obj = Mask()

## filenames of test data
data_filenames = glob.glob(os.path.dirname(pyechomask.__file__).rsplit('\\',1)[0]+'/data/*')

## parse PERG obj and output Sv_dict (see readers.py)
Sv_dict = read_PERGobjs(data_filenames[0],mask_value = mask_obj.mask_value)

## plot 18 kHz echogram
plot_Sv_mask(Sv_dict[18],add_mask = False)

## add masks, defined by method and parameter values
params        = {'fq':18,'th':-75}
threshold_18  = mask_obj.add_mask(binary_threshold,params) ## returns ID of mask
params        = {'fq':18,'mask_value':mask_obj.mask_value}
pulse_18      = mask_obj.add_mask(binary_pulse,params) ## returns ID of mask

## list mask definitions
mask_obj.list_masks()

## EXAMPLE - build pulse/surface mask and apply it (set masked values to mask_value:-999)
Sv_dict[18] = mask_obj.build_mask(Sv_dict,pulse_18,apply = True)

## plot 18 kHz echogram with mask
plot_Sv_mask(Sv_dict[18],add_mask = False)

## EXAMPLE - build single mask (dont apply it)
Sv_dict[18] = mask_obj.build_mask(Sv_dict,threshold_18)

## plot 18 kHz echogram with threshold mask
plot_Sv_mask(Sv_dict[18])

## add more masks
params         = {'fq':38,'th':-85}
threshold_38   = mask_obj.add_mask(binary_threshold,params) ## returns ID of mask
Sv_dict[38]    = mask_obj.build_mask(Sv_dict,threshold_38)
plot_Sv_mask(Sv_dict[38])

## list mask definitions
mask_obj.list_masks()

## EXAMPLE - build composite mask (presence/absence mask)
## e.g. mask if 18 is masked OR 38 is masked
masks             = [threshold_18,threshold_38]
output_mask       = threshold_18
Sv_dict[18]       = mask_obj.build_composite_mask(Sv_dict,masks,output_mask)
plot_Sv_mask(Sv_dict[18])

## EXAMPLE - build composite mask (composite binary mask)
Sv_dict[18] = mask_obj.build_composite_mask(Sv_dict,masks,output_mask,output = 'bitwise')
plot_mask(Sv_dict[18])

'''
In this example, the bitint mask has 4 values (0,1,2,3). 
their binary representations are:
'''
for i in np.unique(Sv_dict[18]['mask']):
    print(i,bin(i)[2:].ljust(2,'0'))    

'''
By example, cells with the value of 3 (11) have values of 1 for the 
first two binary layers.
In this case, the Sv value is larger than -75 dB at 18 and larger
than -85 dB at 38 kHz.
'''
