#%%
import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt

#%%
def get_tglc_phot( filename ):

    # import photometry from TGLC and calculate the magnitudes for PSF and Aperture photometry.
    # Check out: https://github.com/TeHanHunter/TESS_Gaia_Light_Curve/blob/main/tutorial/TGLC_tutorial.ipynb

    """ Example 
 
    filename = '/home/au4148/Dropbox/Undervisning/Master_Projects/Jakob_NGC188/TGLC/V11/TIC461601177/lc/hlsp_tglc_tess_ffi_gaiaid-573941053907094144-s0059-cam3-ccd2_tess_v1_llc.fits'

    import get_tglc_phot

    t, p, a = get_tglc_phot.get_tglc_phot( filename )

    plt.plot( t, a, '.', label = 'Aperture Photometry' )
    plt.plot( t, p, '.', label = 'PSF Photometry' )
    plt.title( 'NGC 188 - V11 - Plotting example' )
    plt.xlabel( 'Time [d]' )
    plt.ylabel( 'Magnitude' )
    plt.legend()

    """




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
#%%
filename =(r'C:\Users\jacob\OneDrive\Dokumenter\Aarhus Universitet\10. Semester\Speciale\Data\TGLC\V11_TGLC\hlsp_tglc_tess_ffi_gaiaid-573941053907094144-s0025-cam4-ccd1_tess_v1_llc.fits')

#import get_tglc_phot

t, p, a = get_tglc_phot( filename )

plt.plot( t, a, '.', label = 'Aperture Photometry' )
plt.plot( t, p, '.', label = 'PSF Photometry' )
plt.title( 'NGC 188 - V11 - Plotting example' )
plt.xlabel( 'Time [d]' )
plt.ylabel( 'Magnitude' )
plt.legend()