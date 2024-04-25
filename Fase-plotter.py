#Fase-plotter
import numpy as np
from astropy.coordinates import SkyCoord
from lightkurve import KeplerTargetPixelFile
import matplotlib.pyplot as plt
from astropy.io import ascii

V18 = 'V18_Rc_data.txt'
Time, Magnitude, MaError = np.loadtxt(V18, unpack=True) #.dat files
#ph, Magnitude, l1, l2, l3, unitmag, rv_a, rv_b = np.loadtxt(V18, unpack=True)


p= 18.798638 #14.469918 # for V20
ph=(Time % p)/p    #for .dat

plt.figure()
plt.clf
plt.plot(ph,-Magnitude,'.')
plt.plot(ph+1,-Magnitude,'.')