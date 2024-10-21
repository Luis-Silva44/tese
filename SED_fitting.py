# %%
from get_flux_values import get_flux_values

#%%
# %% 
def star_values(gaia_id):
    gaia_catalog = "I/355/gaiadr3"  # Gaia DR3 catalog
    gaia_values = Vizier.query_constraints(catalog=gaia_catalog, Source=gaia_id)

    Teff = float(gaia_values[0]['Teff'])
    log_g = float(gaia_values[0]['logg'])
    metalicity = float(gaia_values[0]['__Fe_H_'])

    return Teff, log_g, metalicity

star_values('5707485527450614656')

# %% 
import pysynphot as S
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import LinearNDInterpolator

# %% 

metallicity_grid = np.array([-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.2, 0.5])
logg_grid = np.array([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0])
a = list(range(3000, 13001, 250))
b = list(range(14000, 50001, 1000)) 
Teff_grid = a + b
Teff_grid = np.array(Teff_grid)

def fetch_seds(teff_vals, logg_vals, metallicity_vals, model_name):
    sed_data = []
    
    for teff in teff_vals:
        for logg in logg_vals:
            for metallicity in metallicity_vals:
                try:
                    # Fetch the SED for this grid point
                    sed = S.Icat(model_name, teff, metallicity, logg)
                    sed_data.append(((teff, logg, metallicity), sed))
                except Exception as e:
                    print(f"Error fetching SED for Teff={teff}, logg={logg}, [Fe/H]={metallicity}: {e}")
    
    return sed_data

# %% 
def interpolate_sed(teff, logg, metallicity, model_name):

    # Find the closest lower and upper bounds for interpolation
    teff_low = max([t for t in Teff_grid if t <= teff])
    teff_high = min([t for t in Teff_grid if t > teff])
    logg_low = max([g for g in logg_grid if g <= logg])
    logg_high = min([g for g in logg_grid if g > logg])
    metallicity_low = max([m for m in metallicity_grid if m <= metallicity])
    metallicity_high = min([m for m in metallicity_grid if m > metallicity])
    
    # Fetch SEDs for the surrounding grid points
    sed_data = fetch_seds([teff_low, teff_high], [logg_low, logg_high], [metallicity_low, metallicity_high], model_name)
    
    # Extract wavelengths and fluxes from the fetched SEDs
    wavelengths = sed_data[0][1].wave  # Assuming all SEDs have the same wavelength grid
    fluxes = []
    points = []
    
    for (parameters, sed) in sed_data:
        points.append(parameters)  # (Teff, logg, metallicity)
        fluxes.append(sed.flux)     # Corresponding flux for that point
    
    # Convert to numpy arrays
    points = np.array(points)
    fluxes = np.array(fluxes)
    
    # Create an interpolator for the fluxes at each wavelength
    interpolated_fluxes = []
    for i in range(len(wavelengths)):
        flux_interpolator = LinearNDInterpolator(points, fluxes[:, i])
        interpolated_fluxes.append(flux_interpolator((teff, logg, metallicity)))
   
    interpolated_fluxes = np.array(interpolated_fluxes)
    
    # Create an SED object with the interpolated fluxes
    # interpolated_sed = S.ArraySpectrum(wave=wavelengths, flux=interpolated_fluxes, fluxunits='flam')
    
    return wavelengths, interpolated_fluxes

wavelengths, sed_interpolated = interpolate_sed(4400, 1.5, -0.7, 'ck04models')

plt.plot(wavelengths, sed_interpolated)

