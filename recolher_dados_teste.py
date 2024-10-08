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
    result = Irsa.query_region(coords, catalog='allwise_p3as_psd', spatial='Cone', radius = .3*u.arcmin)
    flux_values = result.to_pandas()
    return flux_values


# %% 
def two_mass_data(star_name):
    coords = Simbad_coords(star_name)
    result = Irsa.query_region(coords, catalog='fp_psc', spatial='Cone', radius = .05*u.arcmin)
    flux_values = result.to_pandas()
    return flux_values


# %% 
def GAIA_data(star_name):
    coords = Simbad_coords(star_name)
    search = Gaia.cone_search_async(coords, radius=u.Quantity(0.7, u.arcmin))
    result = search.get_results()
    flux_values = result.to_pandas()
    return flux_values


# %% 
wise_designation = wise_data("20 LMi")
w1_mag = wise_designation.w1mpro[0]
w2_mag = wise_designation.w2mpro[0]
#w3_mag = wise_designation.w3mpro[0]

# %% 
two_mass_designation = two_mass_data('20 LMi')
j_mag = two_mass_designation.j_m[0]
h_mag = two_mass_designation.h_m[0]
k_mag = two_mass_designation.k_m[0]


print(j_mag)
# %% 
gaia_designation = GAIA_data("20 LMi")
gaia_g_mag = gaia_designation.phot_g_mean_mag[0]
gaia_bp_mag = gaia_designation.phot_bp_mean_mag[0]
gaia_rp_mag = gaia_designation.phot_rp_mean_mag[0]

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
        flux = 2.80e-14 * 10**(-2.5 / mag)
    elif band == 'W2':
        flux = 1.50e-24 * 10**(-2.5 / mag)
    return flux

# %%  
G_flux = flux_from_mag(gaia_g_mag, 'G')
B_flux = flux_from_mag(gaia_bp_mag, 'B')
R_flux = flux_from_mag(gaia_rp_mag, 'R')
J_flux = flux_from_mag(j_mag, 'J')
H_flux = flux_from_mag(h_mag, 'H')
K_flux = flux_from_mag(k_mag, 'K')
W1_flux = flux_from_mag(w1_mag, 'W1')
W2_flux = flux_from_mag(w2_mag, 'W2')

# %% 
flux_values = np.array([G_flux, B_flux, R_flux,J_flux, H_flux, K_flux, W1_flux, W2_flux])
print(flux_values)
plt.plot(flux_values,'o')