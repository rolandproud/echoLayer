#"/usr/bin/env python3
'''

Examples of using binary_impulse mask

Modification History:

'''
## import packages
import matplotlib.pyplot as plt
import pickle
from pyechoplot.plotting import plot_Sv, plot_mask, save_png_plot

## import pyechomask modules
from pyechomask.masks import binary_impulse

## get Sv
filename = './data/krill-Sv38.pkl'
Sv38 = pickle.load( open( filename, 'rb' ) )

## plot the echogram
# plot_Sv(Sv38.transpose())
# plt.title("Krill swarm Sv38")
# plt.show()

## get mask
mask = binary_impulse(Sv38, 10)

# plot_mask(mask.transpose())
# plt.title("Impulse noise mask")
# plt.show()

## plot
plt.figure(1)
plt.subplot(211)
plot_Sv(Sv38.transpose())
plt.subplot(212)
plot_mask(mask.transpose())
plt.show()

# save
#save_png_plot('./','impulse noise')


