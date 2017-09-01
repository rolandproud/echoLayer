# -*- coding: utf-8 -*-
"""
.. :module:: masks
    :synopsis: Contains standardised mask methods
              
              Method format:
 
                def [mask-type]_[unique-name](*args):
                    '''
                    :param [name]: [desc]
                    :type  [name]: [type]
                    
                                   ...
                    return:
                    :param [name]: [desc]
                    :type  [name]: [type]
                    
                    desc: [description]
                    
                    defined by [initials of developer]
                    
                    status: [status(dev,test or product)]
                    
                    '''
                    [code...]
                    return mask
                
             

             mask-type can be 'binary' (0 or 1), 'flag'(range of ints) 
             or 'cont' (continuous: values range from 0-1)
             
             for binary masks: 1 = signal; 0 = noise

| Developed by: Roland Proud (RP) <rp43@st-andrews.ac.uk> 
|               Pelagic Ecology Research Group, University of St Andrews
| Contributors:
|               Rob Blackwell (RB) BAS
| Maintained by:
| Modification History:      
|
UPDATE DESCRIPTIONS AND COMMENTS
"""

## import packages
import numpy as np
from pyechomask.manipulate import median_1D_filter
################################################################## background noise

## background noise removed by readers

################################################################## signal masks 

def binary_threshold(Sv,max_threshold = 999,min_threshold = -999):
    '''
    :param Sv: gridded Sv values (dB re 1m^-1)
    :type  Sv: numpy.array
    
    :param threshold: threshold-value (dB re 1m^-1)
    :type  threshold: float

    return:
    :param mask: binary mask (0 - noise; 1 - signal)
    :type  mask: 2D numpy.array
    
    desc: generate threshold mask
    
    defined by RP
    
    status: product
    
    '''
    ## create mask grid
    mask = np.zeros(Sv.shape).astype(int)
    
    ## apply min threshold
    mask[Sv > min_threshold] = 1
    
    ## apply max threshold
    mask[Sv > max_threshold] = 0
            
    return mask

## detect aggregates/SSLs


def binary_signal(Sv,pl,sample_int,min_sep,max_thickness,max_steps = 10):
    '''
    
    desc: generate signal mask based on:
        Proud R, Cox MJ, Wotherspoon S, Brierley AS. 
        A method for identifying Sound Scattering Layers and extracting key characteristics. 
        Methods Ecol Evol 2015;6:1190–8. doi:10.1111/2041-210X.12396.
    
    defined by RP
    
    status: test
    
    '''

    min_sep       = int(min_sep/sample_int) ## in rows
    max_thickness = int(max_thickness/sample_int) ## in rows
    ## set parameters - limited by processing speed
    ## max size > max thickness of SSL
    ## min size = min seperation distance
    sizes = range(min_sep,int(np.ceil(max_thickness/2 + 1.5*min_sep)),min_sep)
    
    ## linear
    linear  = 10**(Sv/10.)
    linear  = np.ma.masked_invalid(linear)
    
    ## min step distance - set to half of shell length
    step         = int(np.ceil(0.5*pl/1000./(sample_int/1500.))) 
    maxL         = max(sizes)
    
    ## create mirrored image
    row,col = linear.shape
    top     = linear[0:maxL,:][::-1]
    bottom  = linear[row-maxL:,:][::-1]
    image   = np.ma.vstack((top,linear,bottom))
    
    ## create signal mask
    signal = np.zeros(linear.shape)
    
    for size1 in sizes:
        step1      = max([int(np.ceil(size1/max_steps)),step]) 
        sizeRange1 = np.arange(step,int(size1)+1,step1)
        ## build sum of levels above pixel
        mu1 = np.zeros((len(sizeRange1),linear.shape[0],linear.shape[1]))
        for k,l in enumerate(sizeRange1):
            mu1[k,:] = image[maxL-l:row + maxL-l,:]
        for size2 in sizes:
            step2      = max([int(np.ceil(size2/max_steps)),step])  
            sizeRange2 = np.arange(step,int(size2)+1,step2)
            ## build sum of levels below pixel
            mu2 = np.zeros((len(sizeRange2),linear.shape[0],linear.shape[1]))
            for k,l in enumerate(-sizeRange2):
                mu2[k,:] = image[maxL-l:row + maxL-l,:]
            # update signal mask ## MEDIAN!!
            sig = (linear > np.ma.median(mu1,axis = 0)) * \
            (linear > np.ma.median(mu2,axis = 0))  ## where pixel value greater then both upper mean and lower mean then classify as signal
            signal[sig == True] = 1 ## update signal mask
    
    return signal



################################################################### noise masks

## transmit pulse and near-field

def binary_pulse(Sv,noise_level = -999):
    '''
    :param Sv: gridded Sv values (dB re 1m^-1)               
    :type  Sv: numpy.array
    
    :param noise_level: level of background noise (db re 1m^-1)
    :type  noise_level: float

    return:
    :param mask: binary mask (0 - noise; 1 - signal)
    :type  mask: 2D numpy.array
    
    desc: generate pulse mask, mask pulse and surface noise
    
    defined by RP
    
    status: dev
    
    '''
    ## create mask grid
    mask = np.ones(Sv.shape).astype(int)
    
    ## mask pulse and signal up to first noise sample   
    samples,pings = Sv.shape
    for p in range(pings):
        idx           = np.where(Sv[:,p] <= noise_level)[0][0]
        mask[0:idx,p] = 0
         
    return mask


## surface (bubbles/airation)

## seabed

def binary_seabed(Sv, min_depth = 0, threshold = -40, buffer = 5,window_size = 30,noise_level = -999):
    '''
    
    desc: generate seabed mask:

    defined by RP
    
    status: test
    
    '''
    
    ## create mask 
    mask = np.ones(Sv.shape)
    
    ## mask noise
    Sv = np.ma.masked_where(Sv <= noise_level,Sv)
    
    ## mask pulse (if not already removed)
    row,col = Sv.shape
    for c in range(col):
        idx         = np.ma.where(Sv[:,c] < (threshold - 10))[0][0] ## first idx
        Sv[0:idx,c] = noise_level
    
    ## take a guess
    maxidx = np.ma.argmax(Sv,axis = 0)
    
    ## remove values below threshold or above min_depth
    maxidx2 = []
    for k,idx in enumerate(maxidx):
        if Sv[idx,k] < threshold or idx < min_depth:
            maxidx2.append(0)
        else:
            maxidx2.append(idx)
    
    ## if nothing return
    if np.sum(maxidx2) == 0:
        return mask,maxidx2
    
    ## mask bad values        
    maxidx2     = np.array(maxidx2)
    maxidx2     = np.ma.masked_where(maxidx2 == 0,maxidx2)
    
    ## calculate global stats
    #mu = np.ma.mean(maxidx2)
    #sd = np.ma.std(maxidx2)
    ### mask outliers
    #maxidx2 = np.ma.masked_where(maxidx2 > mu + 3*sd,maxidx2)
    #maxidx2 = np.ma.masked_where(maxidx2 < mu - 3*sd,maxidx2)
    
    ## run a local median filter to remove spikes
    maxidx3 = median_1D_filter(maxidx2,window_size)
    
    ## create seabed mask    
    for k,idx in enumerate(maxidx3):
        idx     = int(idx)
        if idx == 0:
            continue
        idx = max([0,idx - buffer])
        mask[idx:,k] = 0
    
    return mask,maxidx3


## false-bottom

## dropped ping/attenuated signal

## transient noise/noise spike

## impulse/interference - regular discrete pulses of sound from external source

def binary_impulse(Sv, threshold):
    '''
    :param Sv: gridded Sv values (dB re 1m^-1)
    :type  Sv: numpy.array
    
    :param threshold: threshold-value (dB re 1m^-1)
    :type  threshold: float
    
    return:
    :param mask: binary mask (0 - noise; 1 - signal)
    :type  mask: 2D numpy.array

    desc: generate threshold mask
    
    defined by RB
    
    status: test
    
    '''
    
    mask = np.ones(Sv.shape).astype(int)

    samples,pings = Sv.shape

    for sample in range(1, samples-1):
        for ping in range(0, pings):
            
            a = Sv[sample-1, ping]
            b = Sv[sample, ping]
            c = Sv[sample+1, ping]

            if (b - a > threshold) & (b - c > threshold):
                mask[sample, ping] = 0

    return mask


## lowered instrument (CTD etc.)







































