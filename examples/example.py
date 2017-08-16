'''

Examples of using masks

Modification History:

'''

## import packages
import os
import glob
import numpy as np

## import user modules
import pyechomask
from pyechomask.masks import binary_pulse, binary_threshold
from pyechomask.readers import get_Sv
from pyechomask.plotting import plot_Sv, plot_mask
from pyechomask.manipulate import merge_binary

## filenames of test data
data_filenames = glob.glob(os.path.dirname(pyechomask.__file__).rsplit('\\',1)[0]+'/data/*')

## parse PERG obj and output Sv (see readers.py)
Sv18 = get_Sv(data_filenames[0],18)
Sv38 = get_Sv(data_filenames[0],38)

## plot 18 kHz echogram
plot_Sv(Sv18)

## create masks
pulse_mask_18     = binary_pulse(Sv18)
threshold_mask_18 = binary_threshold(Sv18,-75)
threshold_mask_38 = binary_threshold(Sv38,-85)

## plot 18 kHz echogram with pulse mask
plot_Sv(Sv18,mask = pulse_mask_18)

#### create composite masks
## presence absence mask
pa_mask              = threshold_mask_18 + threshold_mask_38
pa_mask[pa_mask > 0] = 1
plot_Sv(Sv18,mask = pa_mask)
## merge masks
merged_mask = merge_binary([threshold_mask_18,threshold_mask_38])
## this time, plot just the mask
plot_mask(merged_mask)
#In this example, the merged_mask has 4 values (0,1,2,3). 
#their binary representations are:
for i in np.unique(merged_mask):
    print(i,bin(i)[2:].ljust(2,'0'))    

#By example, cells with the value of 3 (11) have values of 1 for the 
#first two binary mask.
#In this case, the Sv value is larger than -75 dB at 18 and larger
#than -85 dB at 38 kHz.
