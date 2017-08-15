# -*- coding: utf-8 -*-
"""
.. :module:: manipulate
    :synopsis: manipulate masks

| Developed by: Roland Proud (RP) <rp43@st-andrews.ac.uk> 
|               Pelagic Ecology Research Group, University of St Andrews
| Contributors:
|
| Maintained by:
| Modification History:      
|
"""

import numpy as np

def merge_binary(masks):
    '''
    :param masks: list of masks              
    :type  masks: list[numpy.array,...]
    
    :return
    :param output_mask: mask of integers, base2 binary representation of each integer 
                        corresponds to value of each mask input, in the same order.
                        e.g. '0101011', is represented by a unique integer, 43 and 
                        correspond to mask values of 0 for the first mask, 1 for the 
                        second mask etc.
    :type  output_mask: numpy.array of integers
    
    NOTE: the shape of the mask is determined by the first mask in the list
          all masks should have the same number of columns/pings
    
    '''
    ## create array of strings using first mask
    output_mask     = masks[0].astype('S1')
    out_row,out_col = output_mask.shape
    
    for m in masks[1:]:
        row,col                              = m.shape
        mask_new                             = np.zeros((out_row,out_col)) 
        mask_new[0:min(row,out_row),:]       = m[0:min(row,out_row),:]           
        output_mask                          = np.core.defchararray.add\
        (output_mask, mask_new.astype('S1'))
    ## change to integers base 2
    output_mask = np.reshape(np.array([int(x,2) for x in \
                                       output_mask.flatten()]),output_mask.shape)
    
    return output_mask
    
    
    
    