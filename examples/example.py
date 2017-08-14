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

## import user modules
import pyechomask
from pyechomask.mask import Mask
from pyechomask.masks import *
from pyechomask.readers import read_PERGobjs

## filenames of test data
data_filenames = glob.glob(os.path.dirname(pyechomask.__file__).rsplit('\\',1)[0]+'/data/*')

## create instance of mask class
mask_obj = Mask()

## add masks, defined by method and parameter values
params    = {'fq':70,'th':-70}
mask_ID1 = mask_obj.add_mask(binary_threshold,params)
params    = {'fq':120,'th':-70}
mask_ID2 = mask_obj.add_mask(binary_threshold,params)
params    = {'fq':200,'th':-70}
mask_ID3 = mask_obj.add_mask(binary_threshold,params)

## list mask definitions
mask_obj.list_masks()

## EXAMPLE - build single mask
mask = mask_obj.build_mask(read_PERGobjs(data_filenames[0]),mask_ID2)['mask']

## plot mask
plt.figure(figsize = (20,5))
plt.imshow(mask,aspect='auto',cmap = plt.cm.gray) 

## EXAMPLE - build composite mask (presence/absence mask)
## e.g. This would be used to combine all noise masks into a single binary mask
masks       = [mask_ID1,mask_ID2,mask_ID3]
output_mask = mask_ID1
pa_mask     = mask_obj.build_composite_mask(read_PERGobjs(data_filenames[0]),masks,output_mask)

## plot pa mask 
plt.figure(figsize = (20,5))
plt.imshow(pa_mask['mask'],aspect='auto',cmap = plt.cm.gray) 
  
## EXAMPLE - build composite mask (composite binary mask)
bitint_mask = mask_obj.build_composite_mask(read_PERGobjs(data_filenames[0]),masks,output_mask,output = 'bitwise')

## plot bitwise mask 
plt.figure(figsize = (20,5))
plt.imshow(bitint_mask['mask'],aspect='auto',cmap = plt.cm.spectral)
plt.show()

'''
In this example, the bitint_mask has 8 values (0,1,2,3,4,5,6,7). 
their binary representations are:
'''
for i in np.unique(bitint_mask['mask']):
    print(i,bin(i)[2:])    

'''
By example, cells with the value of 6 (110) have values of 1 for the 
first two binary layers and a value of 0 for the last binary layer.
In this case, the Sv value is larger than -70 dB at 70 and 120 kHz but 
smaller than -70 dB at 200 kHz.
'''
