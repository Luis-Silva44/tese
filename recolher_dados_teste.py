# %% 
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.irsa import Irsa
from astroquery.simbad import Simbad
import pandas as pd
from astroquery.gaia import Gaia
import numpy as np
import matplotlib.pyplot as plt

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
def flux_from_mag(mag, band):
    if band == 'G':
        flux = 1.16e-23 * 10**(-2.5 / mag)
    elif band == 'B':
        flux = 4.26e-23 * 10**(-2.5 / mag)
    elif band == 'R':
        flux = 3.08e-23 * 10**(-2.5 / mag)
    elif band == 'J':
        flux = 1.60e-23 * 10**(-2.5 / mag)
    elif band == 'H':
        flux = 1.02e-23 * 10**(-2.5 / mag)
    elif band == 'K':
        flux = 6.40e-24 * 10**(-2.5 / mag)
    elif band == 'W1':
        flux = 2.80e-24 * 10**(-2.5 / mag)
    elif band == 'W2':
        flux = 1.50e-24 * 10**(-2.5 / mag)
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
    wavelen_values =  np.array([230,445,658,1250,1650,2200,3500,4800])
    return flux_values, wavelen_values
# %% 
flux_values, wavelen_values = flux_inital_estimation('Cl* NGC 2632 HSHJ 430')
print(flux_values)
plt.plot(wavelen_values, flux_values,'o')
#%% 
flux_values, wavelen_values = flux_inital_estimation('HD 197481')
print(flux_values)
plt.plot(wavelen_values, flux_values,'o')