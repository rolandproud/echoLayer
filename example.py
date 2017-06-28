'''

Examples of using Mask class and Layer methods


'''
## import os                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
import os

## set wd 
root = 'C:/Users/Roland Proud/Dropbox/py_echolab development/echoLayer/'
os.chdir(root)

## import packages
import pickle    
import matplotlib.pyplot as plt

## import user modules
from Mask import Mask
from layers import *
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
## get example data
with open('testdata', 'rb') as input:
        obj=pickle.load(input)

## create instance of mask class
mask_obj = Mask()

## add layers, defined by method and parameter values
params    = {'fq':38,'th':-80}
layer_ID1 = mask_obj.add_layer(binary_threshold,params)
params    = {'fq':38,'th':-70}
layer_ID2 = mask_obj.add_layer(binary_threshold,params)
params    = {'fq':38,'th':-60}
layer_ID3 = mask_obj.add_layer(binary_threshold,params)

## list layer definitions
mask_obj.list_layers()

## EXAMPLE - build single layer
layer = mask_obj.build_layer(obj,layer_ID1)

## plot layer
plt.figure(figsize = (20,5))
plt.imshow(layer,aspect='auto',cmap = plt.cm.gray)   

## EXAMPLE - build composite layer (presence/absence mask)
## This would be used to combine all noise layers into a single binary layer
layers  = [layer_ID1,layer_ID2,layer_ID3]
pa_mask = mask_obj.build_composite_layer(obj,layers)

## plot pa mask 
plt.figure(figsize = (20,5))
plt.imshow(pa_mask,aspect='auto',cmap = plt.cm.gray) 
  
## EXAMPLE - build composite layer (composite binary mask)
layers      = [layer_ID1,layer_ID2,layer_ID3]
bitint_mask = mask_obj.build_composite_layer(obj,layers,output = 'bitwise')

## plot bitwise mask 
plt.figure(figsize = (20,5))
plt.imshow(bitint_mask,aspect='auto',cmap = plt.cm.jet)

'''
In this example, the bitint_mask has 4 values (0,4,6,7). 
their binary representations are:
'''
for i in np.unique(bitint_mask):
    print(i,bin(i))    

'''
By example, cells with the value of 6 (0b110: '110') have values of 1 for the 
first two binary layers and a value of 0 for the last binary layer.
In this case, the Sv value is larger than -80 and -70 but smaller than -60 dB.


