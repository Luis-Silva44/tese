# %% 
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.irsa import Irsa
from astroquery.simbad import Simbad
import pandas as pd
from astroquery.gaia import Gaia
import numpy as np
import matplotlib.pyplot as plt
import scipy

# %% 
def Simbad_coords(star_name):
    result = Simbad.query_object(star_name)
    star_data = result.to_pandas()

    ra = star_data.RA[0]
    dec = star_data.DEC[0]
    coords = SkyCoord(ra=ra, dec=dec, unit=(u.hourangle, u.deg), frame='icrs')
    return coords
# %% 
def wise_data(star_name):
    coords = Simbad_coords(star_name)
    result = Irsa.query_region(coords, catalog='allwise_p3as_psd', spatial='Cone', radius = 1*u.arcmin)
    flux_values = result.to_pandas()
    return flux_values


# %% 
def two_mass_data(star_name):
    coords = Simbad_coords(star_name)
    result = Irsa.query_region(coords, catalog='fp_psc', spatial='Cone', radius = 1*u.arcmin)
    flux_values = result.to_pandas()
    return flux_values


# %% 
def GAIA_data(star_name):
    coords = Simbad_coords(star_name)
    search = Gaia.cone_search_async(coords, radius=u.Quantity(1, u.arcmin))
    result = search.get_results()
    flux_values = result.to_pandas()
    return flux_values

# %% 
def flux_from_mag(mag, band): # in W / m^2 / micrometers
    if band == 'G':
        flux = 2.82e-8 * 10**(-0.4 * mag)
    elif band == 'B':
        flux = 4.05e-8 * 10**(-0.4 * mag)
    elif band == 'R':
        flux = 1.40e-8 * 10**(-0.4 * mag)
    elif band == 'J':
        flux = 3.129e-13 * 10**(-0.4 * mag)
    elif band == 'H':
        flux = 1.133e-13 * 10**(-0.4 * mag)
    elif band == 'K':
        flux = 4.283e-14 * 10**(-0.4 * mag)
    elif band == 'W1':
        flux = 8.18e-15 * 10**(-0.4 * mag)
    elif band == 'W2':
        flux = 2.42e-15 * 10**(-0.4 * mag)
    return flux

# %% 
def flux_inital_estimation(star_name):
    wise_values = wise_data(star_name)
    w1_mag = wise_values.w1mpro[0]
    w2_mag = wise_values.w2mpro[0]

    two_mass_values = two_mass_data(star_name)
    j_mag = two_mass_values.j_m[0]
    h_mag = two_mass_values.h_m[0]
    k_mag = two_mass_values.k_m[0]

    gaia_values = GAIA_data(star_name)
    gaia_g_mag = gaia_values.phot_g_mean_mag[0]
    gaia_bp_mag = gaia_values.phot_bp_mean_mag[0]
    gaia_rp_mag = gaia_values.phot_rp_mean_mag[0]

    G_flux = flux_from_mag(gaia_g_mag, 'G')
    B_flux = flux_from_mag(gaia_bp_mag, 'B')
    R_flux = flux_from_mag(gaia_rp_mag, 'R')
    J_flux = flux_from_mag(j_mag, 'J')
    H_flux = flux_from_mag(h_mag, 'H')
    K_flux = flux_from_mag(k_mag, 'K')
    W1_flux = flux_from_mag(w1_mag, 'W1')
    W2_flux = flux_from_mag(w2_mag, 'W2')

    flux_values = np.array([G_flux, B_flux, R_flux,J_flux, H_flux, K_flux, W1_flux, W2_flux])
    wavelen_values =  np.array([0.673,0.532,0.797,1.25,1.65,2.15,3.4,4.6]) #in micrometers
    return flux_values, wavelen_values

#%% 
flux_values, wavelen_values = flux_inital_estimation('Vega')
print(flux_values)
plt.plot(wavelen_values, flux_values,'o')

# %% 
from astropy.io import fits

hdul = fits.open('/mnt/c/Users/Luis Silva/Desktop/cenas de faculdade/Mestrado/Tese/ckp00_3500.fits')
hdul.info()
header = hdul[0].header
#print(repr(header))
data = hdul[1].data

# %% 
import pysynphot as ps
import os
os.environ['PYSYN_CDBS'] = '/mnt/c/wsl.localhost/Ubuntu/home/eu/.local/lib/python3.10/site-packages/pysynphot/grid'
sp = ps.Icat('ckp00models', 10000, 0.1, 3.0)
