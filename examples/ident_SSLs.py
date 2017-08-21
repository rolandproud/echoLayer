'''

Use signal filters to identify sound scattering layers
label SSLs and remove small SSLs

Modification History:

'''

## import packages
import matplotlib.pyplot as plt
import gzip
import pickle
from pyechoplot.plotting import plot_Sv, plot_mask, save_png_plot

## import pyechomask modules
from pyechomask.masks import binary_pulse
from pyechomask.manipulate import get_signal_mask, remove_noise, \
                                    signal_row_filter, signal_column_filter, flag

## get Sv data
def getSv(filepath):
    f   = gzip.open(filepath,'rb')
    obj = pickle.load(f,encoding = 'bytes')
    f.close()
    return obj

## set noise level (level of background noise)
## for PERG data, all values below the 95th percentile of the BNL are set to 
## this value
noise_level = -999
    
## read Sv
Sv18 = getSv('./data/PS_Sv18.pklz')

## create signal mask
signal_mask = get_signal_mask(Sv18)

## create pulse mask
pulse_mask_18 = binary_pulse(Sv18)

## remove pulse noise
signal_mask = remove_noise(signal_mask,pulse_mask_18)

## isolate layers
signal_mask = signal_row_filter(signal_mask,25,threshold = 0.8) 
## analysis window set to 25 pings
## therfore, any region of size 1 sample/row by 25 columns/pings that 
## has less than 21 (threshold*window = 20) signal values 
## (values not set to -999 due to low SNR) is set to noise (0), those with more 
## are all set to signal (1)
## window size and threshold dependent upon the SNR and scale of SSLs

## plot 18 kHz echogram with isolated SSLs
plot_Sv(Sv18,mask = signal_mask)
plt.title("18 kHz echogram with pulse removed and sound scattering layers isolated")
plt.show()

## identify pseudo layers - set c. minimun thickness of layers (window = 100 samples)
signal_mask = signal_column_filter(signal_mask,100,threshold = 0.7)

## plot
plt.figure(1)
plt.subplot(211)
plot_Sv(Sv18)
plt.subplot(212)
plot_Sv(Sv18,mask = signal_mask)
plt.title('SSL identification - 18 kHz echosounder data')
plt.show()

## label SSLs and remove small (min_agg_size - number of cells: samples * ping)
labelled_mask = flag(signal_mask,min_agg_size = 5000)

## plot
plt.figure(1)
plt.subplot(211)
plot_Sv(Sv18,mask = signal_mask)
plt.title('flagged SSLs')
plt.subplot(212)
plot_mask(labelled_mask)  
plt.show()

# save
#save_png_plot('./','flagged SSLs')

