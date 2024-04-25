import numpy as np
from astropy.coordinates import SkyCoord
from lightkurve import KeplerTargetPixelFile
import matplotlib.pyplot as plt
from astropy.io import ascii

tpf=KeplerTargetPixelFile("kplr100000941-2009259160929_lpd-targ_2.fits.gz") # read Kepler superstamp target pixel file for Q2 from local disc

# tpf.flux.shape # get format info on targetpixelfile

q=tpf.pipeline_mask # adopt pipeline mask (which is all pixels in aperture for the superstamps)
q[0:len(q)]=False # block all pixels in aperture mask
#q[12:13,77]=True # two-pixel aperture. 
q[17,89]=True # one-pixel aperture. 
#tpf.plot(aperture_mask=q) # plot targetpixelfile with aperture marked on it
lc=tpf.to_lightcurve(aperture_mask=q).flatten() # produce light curve

#lc.plot() # plot flux vs time
#lc.fold(period=18.798638).errorbar() # plot flux vs phase (with period known up front)

#plt.show() # display plots
print(lc.time)
print(lc.flux)

p= 14.469918 #for V20# 18.798638 #for v18#
ph=(lc.time % p)/p

plt.clf
plt.plot(ph,lc.flux,'.')

#ineclipse=[n for n in range(len(lc.flux)) if (ph[n] > 0.330 and ph[n] < 0.355) or (ph[n] > 0.815 and ph[n] < 0.845)]# (ph[n] > 0.330 and ph[n] < 0.355) or (ph[n] > 0.815 and ph[n] < 0.845)]
#outeclipse=[n for n in range(len(lc.flux)) if (ph[n] < 0.330 or ph[n] > 0.845 or (ph[n] > 0.355 and ph[n] < 0.815))] # (ph[n] < 0.330 or ph[n] > 0.845 or (ph[n] > 0.355 and ph[n] < 0.815))

ineclipse=[n for n in range(len(lc.flux)) if (ph[n] > 0.280 and ph[n] < 0.320) or (ph[n] > 0.780 and ph[n] < 0.820)] # (ph[n] > 0.330 and ph[n] < 0.355) or (ph[n] > 0.815 and ph[n] < 0.845)]
outeclipse=[n for n in range(len(lc.flux)) if (ph[n] < 0.280 or ph[n] > 0.820 or (ph[n] > 0.320 and ph[n] < 0.780))] # (ph[n] < 0.330 or ph[n] > 0.845 or (ph[n] > 0.355 and ph[n] < 0.815))
#y=[y[n] for n in range(len(x)) if np.abs(x[n]-RV1) <= (vsini1+10.)]

plt.plot(ph[ineclipse],lc.flux[ineclipse],'.')
plt.plot(ph[outeclipse],lc.flux[outeclipse],'.')
plt.ylabel('Normeret Flux')
plt.xlabel('Fase')

mag=-2.5*np.log10(lc.flux)
rmsmag=np.nanstd(mag[outeclipse])

print(rmsmag)
rmsmag=[rmsmag]*len(mag[ineclipse])
#print(rmsmag)

plt.show()
#%%
plt.clf()
plt.plot(ph[ineclipse],mag[ineclipse],'.')
plt.gca().invert_yaxis()
plt.ylabel('Normeret Flux')
plt.xlabel('Fase')
#plt.xlim(0.3,0.4)

plt.show()
#%%%
A = 'V18_complete_V2.txt'
Time, Magnitude, MaError = np.loadtxt(A, unpack=True) #.dat files

p=  18.798638 #for v18#
ph=(Time % p)/p

plt.clf()
plt.plot(ph,Magnitude,'.')
plt.gca().invert_yaxis()
plt.ylabel('Normeret Flux')
plt.xlabel('Fase')
#plt.xlim(0.3,0.4)
#%%
#Data gennems
time=lc.time[ineclipse]+2454833.0 
mag=mag[ineclipse]                      #Nok ingen grund til det her #%matplotlib qt  extra
rmsmag=np.array(rmsmag)
data = np.array(np.transpose([time,mag,rmsmag]))
ascii.write(data,output='V18_kepler_lc_out_4.dat',names=['BJ-D-TDB','Magnitude','Magnitude error'],overwrite=True)