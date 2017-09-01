'''
Sound Scattering Layer Extraction Method (SSLEM) example

desc: Extract sound scattering layers from echograms 

Reference:
Proud R, Cox MJ, Wotherspoon S, Brierley AS. 
A method for identifying Sound Scattering Layers and extracting key characteristics. 
Methods Ecol Evol 2015;6:1190–8. doi:10.1111/2041-210X.12396.

example by: Roland Proud (RP) <rp43@st-andrews.ac.uk> 
            Pelagic Ecology Research Group, University of St Andrews

Modification History:

'''

## import packages
import matplotlib.pyplot as plt
import gzip
import pickle
from pyechoplot.plotting import plot_Sv, plot_mask, save_png_plot

## import pyechomask modules
from pyechomask.masks import binary_seabed, binary_signal
from pyechomask.manipulate import  signal_row_filter, signal_column_filter, flag,\
        vertical_merge,feature_median, fill_feature_gaps, break_mask, remove_features


## get Sv data
def getSv(filepath):
    f   = gzip.open(filepath,'rb')
    obj = pickle.load(f,encoding = 'bytes')
    f.close()
    return obj

## get Sv values
Sv = getSv('./data/Sv12.pklz')

## plot
plot_Sv(Sv)
plt.title('12kHz echogram from Southern Ocean')
plt.show()

## Sv observation parameters
noise_level   = -999   # background noise level (dB re 1m^-1)
sample_int    = 0.4    # sample interval (m)
pl            = 16.384 # pulse length (ms)
fq            = 12.5   # Frequency

## SSLEM optimization parameters (recomended setting for regional analysis)
min_sep       = 20     # minimum SSL seperation (m):min = pulse length
max_thickness = 300    # Maximun SSL thickness  (m)
min_size      = 100    # Minimum duration of SSL (pings)
min_thickness = 50     # Minimun thickness of SSL (m)

## thresholds 
minSv = -90
maxSv = -50

############################################################################### SSLEM
## NOTE: the order of these steps is important

## calculate number of rows/samples from meters
min_thickness_rows = int(min_thickness/sample_int)
min_sep_rows       = int(min_sep/sample_int)

## remove pulse 
Sv[0:100,:] = noise_level  

## remove weak signal
Sv[Sv < minSv] = noise_level

## get seabed mask, window_size refers to a rolling median filter
## seabed mask in development - check
seabed_mask, seabed_idx = binary_seabed(Sv,buffer = min_sep_rows,window_size = 10)

## plot
plt.figure()
plt.subplot(311)
plot_Sv(Sv)
plt.title('12kHz echogram seabed masked')
plt.subplot(312)
plot_Sv(Sv,mask = seabed_mask)
plt.show()

## apply seabed mask
Sv[seabed_mask == 0] = noise_level

## mask any seabed spikes/edges
Sv[Sv > maxSv] = noise_level

## plot
plt.figure()
plot_Sv(Sv)
plt.title('12kHz echogram seabed masked')
plt.show()

## SSLEM: identify signal pixels
## this step is the bottle neck
signal_mask = binary_signal(Sv,pl,sample_int,min_sep,max_thickness,max_steps = 10)

## plot
plt.figure()
plot_Sv(Sv,mask = signal_mask)
plt.title('12kHz echogram signal')
plt.show()
   
## row smooth signal mask:
## where majority (threshold = 0.5) pixels are signal make all signal, otherwise label all noise
signal = signal_row_filter(signal_mask,min_size,threshold = 0.5) 
   
## remove SSLs that are too thin
signal = signal_column_filter(signal,min_thickness_rows,threshold = 1)

## plot
plt.figure()
plot_Sv(Sv,mask = signal)
plt.title('12kHz echogram filtered')
plt.show()
  
## label SSls (flag) and remove small SSLs (size < min_size*min_thickness_rows)
signal = flag(signal,min_size*min_thickness_rows)

## merge SSLs that are closer together than the min seperation
signal[signal > 0] = 1
signal             = vertical_merge(signal,min_sep_rows)

## fill internal gaps - max size = min_size*min_thickness_rows
signal = fill_feature_gaps(signal,min_size*min_thickness_rows)

## plot
plt.figure()
plot_mask(signal)
plt.title('12kHz echogram SSL mask merged')
plt.show()

## break merged SSLs into individual SSLs
signal = break_mask(signal)

## plot
plt.figure()
plot_mask(signal)
plt.title('12kHz echogram SSL mask')
plt.show()

## remove small SSLs
signal = remove_features(signal.astype(int),min_size*min_thickness_rows)

## get median values of each SSL
Sv_median = feature_median(Sv,signal)

## plot
plt.figure()
plt.subplot(311)
plot_Sv(Sv)
plt.title('12kHz echogram SSLEM: echogram, SSL mask and SSL median')
plt.subplot(312)
plot_mask(signal) 
plt.subplot(313)
plot_Sv(Sv_median)  
plt.show()
save_png_plot('./plots/','SSLEMexample')


