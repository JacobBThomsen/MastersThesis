import numpy as np
from astropy.coordinates import SkyCoord
from lightkurve import KeplerTargetPixelFile
import matplotlib.pyplot as plt

tpf=KeplerTargetPixelFile("kplr100000941-2010265121752_lpd-targ_1.fits.gz") # read Kepler superstamp target pixel file for Q2 from local disc

# tpf.flux.shape # get format info on targetpixelfile

q=tpf.pipeline_mask # adopt pipeline mask (which is all pixels in aperture for the superstamps)
q[0:len(q)]=False # block all pixels in aperture mask
q[2,81]=True # two-pixel aperture. 

tpf.plot(aperture_mask=q) # plot targetpixelfile with aperture marked on it
lc=tpf.to_lightcurve(aperture_mask=q).flatten() # produce light curve

np.savetxt('Targ_1_f.txt', np.transpose([lc.flux])) #Txt af flux file laves
np.savetxt('Targ_1_t.txt', np.transpose([lc.time])) #Txt af flux file laves


lc.plot() # plot flux vs time
lc.fold(period=18.798638).errorbar() # plot flux vs phase (with period known up front)

plt.show() # display plots

#%matplotlib qt          %matplotlib inline

