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
TODO:

UPDATE DESCRIPTIONS AND COMMENTS
"""

import numpy as np
from scipy import ndimage


def median_1D_filter(data,window_size,error_value = 0):
    '''
    Running 1D median filter, width = window_size
    Return error_value where less than 3 data values
    '''
    ## make window odd
    if window_size%2 == 0:
        window_size = window_size +1 ## odd
    ## set min size of window to 3
    window_size = max([window_size,3]) ## min
    ## size of window to the left/right of evaluted pixel
    size   = int(window_size/2)
    ## add mask
    data   = np.ma.masked_invalid(data)  
    result = []
    ## calculate median for rolling window
    ## if  less than 3 data values in window return error_value
    for i in range(size):        
        window      = data[0:size+i]
        window_data = window[window.mask == False]
        if len(window_data) > 2:
            result.append(np.median(window_data))
        else:
            result.append(error_value)
    for i in range(len(data) - 2*size):
        window      = data[i:i+2*size+1]
        window_data = window[window.mask == False]
        if len(window_data) > 2:
            result.append(np.median(window_data))
        else:
            result.append(error_value)           
    for i in range(size):
        window      = data[-2*size+i:]
        window_data = window[window.mask == False]
        if len(window_data) > 2:
            result.append(np.median(window_data))
        else:
            result.append(error_value)
        
    return np.array(result)


def feature_median(Sv,mask,noise_level = -999):
    '''
    for each flagged mask component, calculates median Sv value
    '''
    Sv = np.ma.masked_where(Sv == noise_level,Sv)
    Sv_median = np.ones(mask.shape) * noise_level
    
    for label in np.unique(mask)[1:]:
        idx = np.where(mask == label)
        Sv_median[idx] = np.ma.median(Sv[idx])
        
    return Sv_median


def fill_feature_gaps(mask,max_gap_size = 1000):
    '''
    fill internal gaps of features up to a max size of 
    max_gap_size (in pixels)
    '''
    invert_mask          = np.zeros(mask.shape)
    invert_mask[mask==0] = 1
    invert_mask[mask==1] = 0
    ## fill gaps
    invert_mask = flag(invert_mask,max_gap_size)
    mask[invert_mask == 0] = 1
    
    return mask

def vertical_merge(mask,min_sep):
    '''
    Merge features where distance (in pixels) is less than
    min_sep
    '''
    size    = int(min_sep)
    row,col = mask.shape
    mask2   = np.zeros(mask.shape)
    for asize in range(1,size):
        bsize = size - asize
        above = np.zeros((row-2*size,col))
        below = np.zeros((row-2*size,col))
        for i in range(asize):
            above += mask[size-i-1:row-size-i-1,:]
        for i in range(bsize):   
            below += mask[size + 1 + i:row-size+i+1,:]            
        above[above > 0] = 1 
        below[below > 0] = 1
        mask2[size:row-size,:][(above + below) == 2] = 1
    
    return mask2


def label_ping(ping_mask):
    '''
    label (flag) each seperate feature, 
    defined as signal regions divided by noise, of a single ping
    '''
    ## blank
    ping_label = np.zeros(ping_mask.shape)
    
    ## signal
    idx = np.where(ping_mask > 0)[0]
    
    ## if no signal exit
    if len(idx) == 0:
        return ping_label
    
    ## label each feature
    label = 1
    vals  = [idx[0]]
    for k,i in enumerate(idx[1:]):
        if i != idx[k] + 1:
            ping_label[vals] = label
            label += 1
            vals = [i]
        else:
            vals.append(i)
    ## last one
    ping_label[vals] = label
    return ping_label


def break_mask(mask):
    '''
    break a mask into individual features (no vertical gaps in signal)
    '''
    
    mask[mask > 0]     = 1
    row,col            = mask.shape
    prev_col           = label_ping(mask[:,0])
    labelled_mask      = np.zeros(mask.shape)
    labelled_mask[:,0] = prev_col 
    for c in range(1,col):
        next_col        = label_ping(mask[:,c])
        all_connections = []
        for next_feature in np.unique(next_col)[1:]:
            connections = []
            next_idx = np.where(next_col == next_feature)[0]
            for prev_feature in np.unique(prev_col)[1:]:
                prev_idx = np.where(prev_col == prev_feature)[0]
                for i in next_idx:
                    if i in prev_idx:
                        connections.append(prev_feature)
                        break
            all_connections.append(connections)
        for k,con in enumerate(all_connections):
            new = False
            if len(con) == 1:
                for k2,con2 in enumerate(all_connections):
                    if k == k2:
                        continue
                    if con[0] in con2:
                        new = True
            else:## new connection
                new = True
            idx = np.where(next_col == k+1)[0]
            if new:            
                labelled_mask[idx,c] = np.max(labelled_mask) + 1
            else:
                labelled_mask[idx,c] = con[0]
        prev_col = labelled_mask[:,c]
        
    return labelled_mask


def flag(mask,min_agg_size = 0,struct = None):
    """
    remove small aggregates and label others
    """
    structure=[[1,1,1],
               [1,1,1],
               [1,1,1]]
    if struct != None:
        structure == struct
    
    ## label image
    label_im, nb_labels    = ndimage.label(np.asfarray(mask),structure)
    
    return remove_features(label_im, min_agg_size)

def remove_features(label_im, min_agg_size = 0):
    '''
    remove masked features smaller than min_agg_size (in pixels)
    '''
    mask                   = np.zeros(label_im.shape)
    mask[label_im > 0]     = 1
    sizes                  = ndimage.sum(mask, label_im, range(len(np.unique(label_im))))
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

def signal_column_filter(mask,window,threshold = 0.5):
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
    mask_sum[mask_sum < (window*threshold)] = 0  # noise
    mask_sum[mask_sum >= (window*threshold)]  = 1  # signal
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
    
    
    
    