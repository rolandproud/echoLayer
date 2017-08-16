# -*- coding: utf-8 -*-
"""
.. :module:: plotting
    :synopsis: plotting functions

| Developed by: Roland Proud (RP) <rp43@st-andrews.ac.uk> 
|               Pelagic Ecology Research Group, University of St Andrews
| Contributors:
|
| Maintained by:
| Modification History:      
|
"""
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt


def plot_mask(mask):
    '''
    '''
    ## shape
    row,col = mask.shape
    
    ## plot
    #f, (ax1) = plt.subplots(1, figsize = (20,10))
    p1       = plt.imshow(mask, cmap = plt.cm.spectral,\
            interpolation='nearest',aspect='auto')
    plt.colorbar(p1,pad = 0)
    plt.xlabel('columns',fontsize = 18)
    plt.ylabel('rows',fontsize = 18)

def plot_Sv(Sv,mask = None):
    '''
    '''
    ## get echoview colormap
    setup_ek500_cmap()
    ek500_cmap = mpl.cm.get_cmap('ek500')
    ek500_norm = mpl.colors.BoundaryNorm(np.linspace(-89,-34,12), 12, clip=False)
    
    ## add mask
    if mask is not None:
        Sv = np.ma.masked_where(mask == 0,Sv)
    
    ## shape
    row,col = Sv.shape
    
    ## plot
    #f, (ax1) = plt.subplots(1, figsize = (20,10))
    p1       = plt.imshow(Sv, cmap = ek500_cmap,norm = ek500_norm,\
            interpolation='nearest',aspect='auto')
    plt.colorbar(p1,pad = 0)
    plt.xlabel('columns',fontsize = 18)
    plt.ylabel('rows',fontsize = 18)

def setup_ek500_cmap():
    '''
    Creates the ek500 colormap and boundary norm
    Taken from pyecholab (Add ref)
    '''
    float_ek500_cmap_colorlist = [(0.62, 0.62, 0.62),
                            (0.37, 0.37, 0.37),
                            (0.0, 0.0, 1.0),
                            (0.0, 0.0, 0.498),
                            (0.0, 0.749, 0.0),
                            (0.0, 0.498, 0.0),
                            (1.0, 1.0, 0.0),
                            (1.0, 0.498, 0.0),
                            (1.0, 0.0, 0.749),
                            (1.0, 0.0, 0.0),
                            (0.651, 0.325, 0.235),
                            (0.471, 0.235, 0.157)]

    ek500_cmap = \
        mpl.colors.ListedColormap(float_ek500_cmap_colorlist, name='ek500')
    ek500_cmap.set_bad(color='k', alpha=1.0)
    ek500_cmap.set_under('w', alpha=1.0)
    ek500_cmap.set_over(color=float_ek500_cmap_colorlist[-1], alpha=1.0)
    
    if 'ek500' not in mpl.cm.cmap_d:
        mpl.cm.register_cmap('ek500', cmap=ek500_cmap)
    else:
        mpl.cm.cmap_d['ek500'] = ek500_cmap
