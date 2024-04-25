#%%
import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt
#%% Data
a18 = fits.getdata( 'pathos_tess_lightcurve_tic-00461618602-s018_tess_v1_llc.fits' )
a20 = fits.getdata( 'pathos_tess_lightcurve_tic-00461618602-s020_tess_v1_llc.fits' )
a25 = fits.getdata( 'pathos_tess_lightcurve_tic-00461618602-s025_tess_v1_llc.fits' )
a26 = fits.getdata( 'pathos_tess_lightcurve_tic-00461618602-s026_tess_v1_llc.fits' )
a40 = fits.getdata( 'pathos_tess_lightcurve_tic-00461618602-s040_tess_v1_llc.fits' )
a52 = fits.getdata( 'pathos_tess_lightcurve_tic-00461618602-s052_tess_v1_llc.fits' )
a53 = fits.getdata( 'pathos_tess_lightcurve_tic-00461618602-s053_tess_v1_llc.fits' )
#%% Data extract
per = 6.5043035 #period 

ph = 'best_phot_flux_raw'
ph = 'best_phot_flux_cor'

m18 = 25.0 - 2.5*np.log10( a18[ ph ] )
m20 = 25.0 - 2.5*np.log10( a20[ ph ] )
m25 = 25.0 - 2.5*np.log10( a25[ ph ] )
m26 = 25.0 - 2.5*np.log10( a26[ ph ] )
m40 = 25.0 - 2.5*np.log10( a40[ ph ] )
m52 = 25.0 - 2.5*np.log10( a52[ ph ] )
m53 = 25.0 - 2.5*np.log10( a53[ ph ] )

ph18 = np.mod( a18['time'], per )
ph20 = np.mod( a20['time'], per )
ph25 = np.mod( a25['time'], per )
ph26 = np.mod( a26['time'], per )
ph40 = np.mod( a40['time'], per )
ph52 = np.mod( a52['time'], per )
ph53 = np.mod( a53['time'], per )

i18 = np.where( (ph18 > 4.1 ) & (ph18 < 4.4) )[0]
i20 = np.where( (ph20 > 4.1 ) & (ph20 < 4.4) )[0]
i25 = np.where( (ph25 > 4.1 ) & (ph25 < 4.4) )[0]
i26 = np.where( (ph26 > 4.1 ) & (ph26 < 4.4) )[0]
i40 = np.where( (ph40 > 4.1 ) & (ph40 < 4.4) )[0]
i52 = np.where( (ph52 > 4.1 ) & (ph52 < 4.4) )[0]
i53 = np.where( (ph53 > 4.1 ) & (ph53 < 4.4) )[0]
#%%Plot it all %matplotlib qt and %matplotlib inline
plt.plot( ph18, m18 - np.mean(m18[i18]), '.k',  label = 'Sector 18' )
plt.plot( ph20, m20 - np.mean(m20[i20]), '.y' , label = 'Sector 20')
plt.plot( ph25, m25 - np.mean(m25[i25]), '.m' , label = 'Sector 25')
plt.plot( ph26, m26 - np.mean(m26[i26]), '.b' , label = 'Sector 26')
plt.plot( ph40, m40 - np.mean(m40[i40]), '.c' , label = 'Sector 40')
plt.plot( ph52, m52 - np.mean(m52[i52]), '.g' , label = 'Sector 52')
plt.plot( ph53, m53 - np.mean(m53[i53]), '.r' , label = 'Sector 53')

plt.ylim( 0.5, -0.5 )
plt.grid( ls = 'dotted' )
plt.legend()
#%% Plot some
plt.plot( ph18, m18 - np.mean(m18[i18]), '.k',  label = 'Sector 18' )

plt.ylim( 0.5, -0.5 )
plt.grid( ls = 'dotted' )
plt.legend()