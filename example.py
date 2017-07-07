'''

Examples of using Mask class and Layer methods

Modification History

7/7/2017: added multifrequency example
          added vector representation of layer
'''

## import os
import os

## set wd 
root = 'C:/Users/Roland Proud/Dropbox/py_echolab development/echoLayer/'
os.chdir(root)

## import packages
import pickle    
import matplotlib.pyplot as plt
#import xmltodict
#from xml.etree import ElementTree

## import user modules
from Mask import Mask
from layers import *

## change wd to data directory
os.chdir(root + 'test data/')
       
## read raw multi-frequency EK60 data
with open('rawData', 'rb') as input:
        raw = pickle.load(input)
                                                                                      
## create instance of mask class
mask_obj = Mask()

## add layers, defined by method and parameter values
params    = {'fq':'38','th':-70}
layer_ID1 = mask_obj.add_layer(binary_threshold,params)
params    = {'fq':'120','th':-70}
layer_ID2 = mask_obj.add_layer(binary_threshold,params)
params    = {'fq':'200','th':-70}
layer_ID3 = mask_obj.add_layer(binary_threshold,params)

## list layer definitions
mask_obj.list_layers()

## EXAMPLE - build single layer
layer = mask_obj.build_layer(raw,layer_ID2)

## plot layer
plt.figure(figsize = (20,5))
plt.imshow(layer,aspect='auto',cmap = plt.cm.gray) 

## EXAMPLE - get vector representation of single layer
## less memory to store and will be able to merge and sum layers
## of any size (grid spacing)
## ptv: ping-time-vector, value: [ping range]: [sample range (seconds)]
ptv  = mask_obj.build_layer(raw,layer_ID3,output = 'ptv')

## update build_composite function to read ptv dictionaries


## EXAMPLE - build composite layer (presence/absence mask)
## Currently only works with layers with same grid size...
## e.g. This would be used to combine all noise layers into a single binary layer
layers  = [layer_ID1,layer_ID2,layer_ID3]
pa_mask = mask_obj.build_composite_layer(raw,layers)

## plot pa mask 
plt.figure(figsize = (20,5))
plt.imshow(pa_mask,aspect='auto',cmap = plt.cm.gray) 
  
## EXAMPLE - build composite layer (composite binary mask)
layers      = [layer_ID1,layer_ID2,layer_ID3]
bitint_mask = mask_obj.build_composite_layer(raw,layers,output = 'bitwise')

## plot bitwise mask 
plt.figure(figsize = (20,5))
plt.imshow(bitint_mask,aspect='auto',cmap = plt.cm.spectral)
plt.show()

'''
In this example, the bitint_mask has 8 values (0,1,2,3,4,5,6,7). 
their binary representations are:
'''
for i in np.unique(bitint_mask):
    print(i,bin(i)[2:])    

'''
By example, cells with the value of 6 (110) have values of 1 for the 
first two binary layers and a value of 0 for the last binary layer.
In this case, the Sv value is larger than -70 dB at 38 and 120 kHz but 
smaller than -70 dB at 200 kHz.
'''
