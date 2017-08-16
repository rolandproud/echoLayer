#"/usr/bin/env python3

import matplotlib.pyplot as plt
import pickle

from pyechomask.masks import mask_impulse_noise
from pyechomask.plotting import plot_Sv, plot_mask

filename = '../data/krill-Sv38.pkl'
Sv38 = pickle.load( open( filename, 'rb' ) )

## plot the echogram
# plot_Sv(Sv38.transpose())
# plt.title("Krill swarm Sv38")
# plt.show()

mask = mask_impulse_noise(Sv38, 10)

# plot_mask(mask.transpose())
# plt.title("Impulse noise mask")
# plt.show()


plt.figure(1)
plt.subplot(211)
plot_Sv(Sv38.transpose())
plt.subplot(212)
plot_mask(mask.transpose())
plt.show()
