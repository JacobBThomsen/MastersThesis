import numpy as np
from astropy.coordinates import SkyCoord
from lightkurve import KeplerTargetPixelFile
import matplotlib.pyplot as plt

tpf=KeplerTargetPixelFile("/mnt/hgfs/vmware-shared/NGC6791/V56/Kepler/kplr100000941-2009166043257_lpd-targ.fits.gz") # read Kepler superstamp target pixel file for Q2 from local disc

# tpf.flux.shape # get format info on targetpixelfile

q=tpf.pipeline_mask # adopt pipeline mask (which is all pixels in aperture for the superstamps)
q[0:len(q)]=False # block all pixels in aperture mask
#q[12:13,77]=True # two-pixel aperture. 
q[1,77]=True # one-pixel aperture. 
#tpf.plot(aperture_mask=q) # plot targetpixelfile with aperture marked on it
lc=tpf.to_lightcurve(aperture_mask=q).flatten() # produce light curve

#lc.plot() # plot flux vs time
#lc.fold(period=18.798638).errorbar() # plot flux vs phase (with period known up front)

#plt.show() # display plots
print(lc.time)
print(lc.flux)

p=18.798638
ph=(lc.time % p)/p

plt.clf
plt.plot(ph,lc.flux,'.')

ineclipse=[n for n in range(len(lc.flux)) if (ph[n] > 0.330 and ph[n] < 0.355) or (ph[n] > 0.815 and ph[n] < 0.845)]
outeclipse=[n for n in range(len(lc.flux)) if (ph[n] < 0.330 or ph[n] > 0.845 or (ph[n] > 0.355 and ph[n] < 0.815))]
#y=[y[n] for n in range(len(x)) if np.abs(x[n]-RV1) <= (vsini1+10.)]

plt.plot(ph[ineclipse],lc.flux[ineclipse],'.')
plt.plot(ph[outeclipse],lc.flux[outeclipse],'.')

mag=-2.5*np.log10(lc.flux)
rmsmag=np.std(mag[outeclipse])
print(rmsmag)
rmsmag=[rmsmag]*len(mag)
print(rmsmag)

plt.show()
plt.clf()
plt.plot(ph[ineclipse],mag[ineclipse],'o')

plt.show()

