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
from scipy import ndimage

def flag(mask,min_agg_size = 0):
    """
    remove small aggregates and label others
    """

    structure=[[1,1,1],
              [1,1,1],
              [1,1,1]]

    label_im, nb_labels    = ndimage.label(np.asfarray(mask),structure)
    sizes                  = ndimage.sum(mask, label_im, range(nb_labels + 1))
    mask_size              = sizes < min_agg_size
    remove_pixel           = mask_size[label_im]
    label_im[remove_pixel] = 0
    labels                 = np.unique(label_im)
    label_clean            = np.searchsorted(labels, label_im)

    return label_clean


def signal_row_filter(mask,window,threshold = 0.5):
    '''
    :param Sv: gridded Sv values (dB re 1m^-1)
    :type  Sv: numpy.array
    
    :param window: size in numbwer of pings/columns of analysis window
    :type  window: int
    
    :param threshold: value between 0 and 1, proportion of window required 
                        to assign as signal value
    :type  threshold: float
    
    :param noise_level: level of background noise (db re 1m^-1)
    :type  noise_level: float

    desc: isolate signal rows
    
    defined by RP
    
    status: dev
    '''
    window        = int(window)
    row,col       = mask.shape
    filtered_mask = np.zeros((mask.shape))
    mask_sum      = np.zeros((row,col-window+1))
    for i in range(window):
        mask_sum += mask[:,i:col-(window-i-1)]
    ##
    mask_sum[mask_sum <= (window*threshold)] = 0  # noise
    mask_sum[mask_sum > (window*threshold)]  = 1  # signal
    for i in range(window):
        filtered_mask[:,i:col-(window-i-1)] += mask_sum
    filtered_mask[filtered_mask > 0] = 1
        
    #for r in range(row):
    #    for c in range(int(window/2),int(col-window/2) +1,1):
    #        col_weight = sum(mask[r,c-int(window/2):c+int(window/2)+1])
    #        if col_weight > (window*threshold):
    #            filtered_mask[r,c-int(window/2):c+int(window/2)+1] = 1
                
    return filtered_mask

def signal_column_filter(mask,window,threshold = 0.5,erosion = True):
    '''
    :param Sv: gridded Sv values (dB re 1m^-1)
    :type  Sv: numpy.array
    
    :param window: size in numbwer of samples/rows of analysis window
    :type  window: int
    
    :param threshold: value between 0 and 1, proportion of window required 
                        to assign as signal value
    :type  threshold: float
    
    :param noise_level: level of background noise (db re 1m^-1)
    :type  noise_level: float

    desc: isolate signal rows
    
    defined by RP
    
    status: dev
    '''
    #if erosion:
    window        = int(window) 
    row,col       = mask.shape
    filtered_mask = np.zeros((mask.shape))
    mask_sum      = np.zeros((row-window+1,col))
    for i in range(window):
        mask_sum += mask[i:row-(window-i-1),:]
    ## 
    mask_sum[mask_sum <= (window*threshold)] = 0  # noise
    mask_sum[mask_sum > (window*threshold)]  = 1  # signal
    for i in range(window):
        filtered_mask[i:row-(window-i-1),:] += mask_sum
    filtered_mask[filtered_mask > 0] = 1
    ## extend
    #top      = mask_sum[0:int(window/2),:]
    #bottom   = mask_sum[(row - window - int(window/2)):,:]
    #mask_sum = np.vstack((top,mask_sum,bottom))

    
    ## create mask
    #filtered_mask[mask_sum > (window*threshold)] = 1
    #else:
    #    row,col       = mask.shape
    #    filtered_mask = np.zeros((mask.shape))
    #    for c in range(col):
    #        for r in range(int(window/2),int(row-window/2) +1,1):
    #            row_weight = sum(mask[r-int(window/2):r+int(window/2)+1,c])
    #            if row_weight > (window*threshold):
    #                filtered_mask[r-int(window/2):r+int(window/2)+1,c] = 1
                
    return filtered_mask


def remove_noise(mask,noise_mask):
    '''
    '''
    mask[noise_mask == 0] = 0
    return mask

def get_signal_mask(Sv,noise_level = -999):
    '''
    '''
    signal_mask                    = np.ones((Sv.shape))
    signal_mask[Sv == noise_level] = 0
    return signal_mask


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
    
    
    
    