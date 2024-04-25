
# conda activate tglc
# Check out: https://github.com/TeHanHunter/TESS_Gaia_Light_Curve/blob/main/tutorial/TGLC_tutorial.ipynb

from astroquery.mast import Observations

#ticid  = '461618602'     # V12 - NGC 188 
#per    = 6.5042969      
#target = 'TIC461618602'  # V12 - NGC 188 
 
ticid  = '461601177'      # V11 - NGC 188 :   Ra, Dec = 00 45 22.74  +85 12 38.7)
per    = 35.17756
target = 'TIC461601177'   # V11 - NGC 188 :   Ra, Dec = 00 45 22.74  +85 12 38.7)



import os
import numpy as np
import matplotlib.pyplot as plt
from tglc.quick_lc import tglc_lc


local_directory = f'{target}/'    # directory to save all files

os.makedirs(local_directory, exist_ok=True)

"""
For V11 the following sectors are available (27/2-2023): 18,19,20,25,26,40,52,53,59 
"""

tglc_lc(target=target, 
        local_directory=local_directory, 
        size=90,                # FFI cutsize. Recommand at least 50 or larger for better performance. Cannot exceed 99. 
                                # Downloading FFI might take longer (or even cause timeouterror) for larger sizes. 
        save_aper=False,        # whether to save 5*5 pixels timeseries of the decontaminated images in fits file primary HDU
        limit_mag=15,           # the TESS magnitude lower limit of stars to output
        get_all_lc=False,       # whether to return all lcs in the region. If False, return the nearest star to the target coordinate
        first_sector_only=True, # whether to return only lcs from the sector this target was first observed. 
                                # If False, return all sectors of the target, but too many sectors could be slow to download.
        sector=None,            # If first_sector_only=True, sector will be ignored.
                                # If first_sector_only=False and sector = None, return all observed sectors
                                # If first_sector_only=False and sector = 1, return only sector 1. 
                                # (Make sure only put observed sectors. All available sectors are printed in the sector table.)
        prior=None)             # If None, does not allow all field stars to float. SUGGESTED for first use. 
                                # If float (usually <1), allow field stars to float with a Gaussian prior with the mean 
                                # at the Gaia predicted value the width of the prior value multiplied on the Gaia predicted value.)


##



def get_tglc_phot( filename ):

    from astropy.io import fits 

    hdul         = fits.open( filename )
    hdul.info()

    filt              = list(hdul[1].data['TESS_flags'] == 0) and list(hdul[1].data['TGLC_flags'] == 0)
    # filter out bad datapoints from both TESS FFI flags and TGLC flags

    time           = hdul[1].data['time'][filt]          # Time ... BJD ??

    psf_flux       = hdul[1].data['psf_flux'][filt]      # raw psf flux
    psf_mag        = 25.0 - 2.5*np.log10( psf_flux )
    psf_flux_err   = hdul[1].header['PSF_ERR']           # raw psf flux error
    psf_mag_err    = -1.0857*psf_flux_err/psf_flux       # raw aper flux error

    aper_flux      = hdul[1].data['aperture_flux'][filt] # raw aper flux
    aper_mag       = 25.0 - 2.5*np.log10( aper_flux )
    aper_flux_err  = hdul[1].header['APER_ERR']          # raw aper flux error
    aper_mag_err   = -1.0857*aper_flux_err/aper_flux     # raw aper flux error

    return time, psf_mag, aper_mag 


    #plt.errorbar(time, psf_flux,  psf_flux_err,  '.', label = 'PSF')
    #plt.errorbar(time, aper_flux, aper_flux_err, '.', label = 'Aperture')


    #plt.plot(time, psf_flux,    '.', label = 'PSF')
    #plt.plot(time, aper_flux,  '.', label = 'Aperture')
    #plt.title( target + '- TESS Sector ??' )
    #plt.xlabel('TBJD')
    #plt.ylabel('Flux e-/s')
    #plt.legend()
    #plt.show()
