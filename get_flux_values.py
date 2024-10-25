#%% 
import astropy.units as u
from astroquery.vizier import Vizier
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pysynphot as S
import os

#STILL NEEDED: CORRECT ERRORS; REMOVE WARNING FROM COORDINATES
# %% 

def vizier_coords(gaia_id):
    Vizier.ROW_LIMIT = 1  
    gaia_catalog = "I/355/gaiadr3"  # Gaia DR3 catalog

    gaia_values = Vizier.query_constraints(catalog=gaia_catalog, Source=gaia_id)
    
    if gaia_values:
        ra, dec = gaia_values[0]['RA_ICRS'][0], gaia_values[0]['DE_ICRS'][0]
        return ra, dec
    else:
        raise ValueError("No Gaia id found.")

#%% 
def wise_values(gaia_id):
    wise_catalog = 'II/311/wise'
    ra, dec = vizier_coords(gaia_id)
    wise_values = Vizier.query_region(f"{ra} {dec}", radius=10* u.arcsec , catalog=wise_catalog)
    if wise_values:
        return wise_values
    else:
        raise ValueError('No WISE data found')

def two_mass_values(gaia_id):
    two_mass_catalog = 'II/246/out'
    ra, dec = vizier_coords(gaia_id)

    two_mass_values = Vizier.query_region(f"{ra} {dec}", radius=10*u.arcsec, catalog=two_mass_catalog)
    
    if two_mass_values:
        return two_mass_values
    else: 
        raise ValueError('No 2MASS data found')

def gaia_values(gaia_id):
    gaia_catalog = "I/355/gaiadr3"  # Gaia DR3 catalog
    gaia_values = Vizier.query_constraints(catalog=gaia_catalog, Source=gaia_id)

    return gaia_values

def mag_to_flux(mag,band): # in W cm-2 micrometer-1
    flux_zero_points =  {'J':3.1293e-13,
                         'H':1.133e-13,
                         'K':4.283e-14,
                         'W1':8.180e-15,
                         'W2':2.415e-15
                         }
    flux_constant = flux_zero_points.get(band)

    if flux_constant is None:
        raise ValueError(f'No zero point flux value found for {band} band')
    return flux_constant * 10**(-0.4 * mag)

#%%
def get_flux_values(gaia_id):
    wise_data = wise_values(gaia_id)
    two_mass_data = two_mass_values(gaia_id)
    gaia_data = gaia_values(gaia_id)
    W1_mag, W1_mag_err = float(wise_data[0]['W1mag']), float(wise_data[0]['e_W1mag'])
    W2_mag, W2_mag_err = float(wise_data[0]['W2mag']), float(wise_data[0]['e_W2mag'])

    J_mag, J_mag_err = float(two_mass_data[0]['Jmag']), float(two_mass_data[0]['e_Jmag'])
    H_mag, H_mag_err = float(two_mass_data[0]['Hmag']), float(two_mass_data[0]['e_Hmag'])
    K_mag, K_mag_err = float(two_mass_data[0]['Kmag']), float(two_mass_data[0]['e_Kmag'])

    W1_flux, W1_flux_err = mag_to_flux(W1_mag,'W1'), mag_to_flux(W1_mag_err,'W1')
    W2_flux, W2_flux_err = mag_to_flux(W2_mag,'W2'), mag_to_flux(W2_mag_err,'W2')
    J_flux, J_flux_err = mag_to_flux(J_mag,'J'), mag_to_flux(J_mag_err,'J')
    H_flux, H_flux_err = mag_to_flux(H_mag,'H'), mag_to_flux(H_mag_err,'H')
    K_flux, K_flux_err = mag_to_flux(K_mag,'K'), mag_to_flux(K_mag_err,'K')

    G_flux, G_flux_err = float(gaia_data[0]['FG']), float(gaia_data[0]['e_FG'])
    GBP_flux, GBP_flux_err = float(gaia_data[0]['FBP']), float(gaia_data[0]['e_FBP'])
    GRP_flux, GRP_flux_err = float(gaia_data[0]['FRP']), float(gaia_data[0]['e_FRP'])

    G_flux, G_flux_err = G_flux * 1.346109e-22, G_flux_err * 1.346109e-22
    GBP_flux, GBP_flux_err = GBP_flux * 3.009167E-22, GBP_flux_err * 3.009167E-22
    GRP_flux, GRP_flux_err = GRP_flux * 1.638483E-22, GRP_flux_err * 1.638483E-22

    flux_values = np.array([G_flux,GBP_flux,GRP_flux, W1_flux,W2_flux,J_flux,H_flux,K_flux])
    flux_err = np.array([G_flux_err,GBP_flux_err,GRP_flux_err,W1_flux_err,W2_flux_err,J_flux_err,H_flux_err,K_flux_err])
    
    
    filter_bands = {'G': 0.673,'GBP':0.532, 'GRP':0.797, 
                        'W1':3.4, 'W2':4.6, 
                        'J':1.25, 'H':1.65, 'K':2.15,}

    filter_wavelen = np.array([d for d in filter_bands.values()])

    return filter_wavelen, flux_values, flux_err

# %% 
wavelen, flux_values, flux_error = get_flux_values('2135550755683407232')

plt.plot(wavelen, flux_values,'o')
plt.xlabel('Wavelenght (μm)')
plt.ylabel('Flux (W / cm-2 μm-1)')
plt.title('Flux values of the star for each filter')

print(flux_values)
print(flux_error)