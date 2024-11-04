# %%
from get_flux_values import *
from scipy.interpolate import LinearNDInterpolator
import pysynphot as S
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from uncertainties import ufloat
# %% 
first_temp_steps = list(range(3000,13001,250))
second_temp_steps = list(range(13000, 50001, 1000)) 
Teff_grid = (first_temp_steps + second_temp_steps) 
logg_grid = [0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0] 
mettalicity_grid = [-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5] 

def create_SEDs(Teff_vals, mettalicity_vals, logg_vals):
    SED_data = []

    for teff in Teff_vals:
        for mett in mettalicity_vals:
            for logg in logg_vals:
                try:
                    sed_values = S.Icat('ck04models',teff,mett,logg)
                    SED_data.append(((teff,mett,logg),sed_values))
                except Exception as e:
                    print(f"Error getting SED values for Teff={teff}, log_g={logg} and mettalicity={mett}")

<<<<<<< HEAD
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

=======
    return SED_data

def SED_high_and_low(Teff,mettalicity,logg):
    Teff_low = max([t for t in Teff_grid if t <= Teff])
    Teff_high = min([t for t in Teff_grid if t > Teff])
    mettalicity_low = max([m for m in mettalicity_grid if m <= mettalicity])
    mettalicity_high = min([m for m in mettalicity_grid if m > mettalicity])
    logg_low = max([l for l in logg_grid if l <= logg])
    logg_high = min([l for l in logg_grid if l > logg])

    Teff_values = [Teff_low, Teff_high]
    mettalicity_values =  [mettalicity_low, mettalicity_high]
    logg_values = [logg_low, logg_high]

    SED_data = create_SEDs(Teff_values,mettalicity_values,logg_values)

    return SED_data

def SED_interpolator(Teff,mettalicity,logg):
    SED_data = SED_high_and_low(Teff,mettalicity,logg)
    SED_wavelen = SED_data[0][1].wave * u.angstrom

    fluxes = []
    points = []

    for (parameters, SED_values) in SED_data: 
        fluxes.append(SED_values.flux)
        points.append(parameters)
    
    fluxes = np.array(fluxes)
    points = np.array(points)
    
    interpolated_fluxes = []

    for i in range(len(SED_wavelen)):
        flux_interpolator = LinearNDInterpolator(points, fluxes[:,i])
        interpolated_fluxes.append(flux_interpolator(Teff,mettalicity,logg))

    interpolated_fluxes = np.array(interpolated_fluxes) * u.erg / u.cm**2 / u.s / u.angstrom

    SED_wavelen = SED_wavelen.to(u.um)
    model_flux_Jy = interpolated_fluxes.to(u.Jy, u.spectral_density(SED_wavelen))

    return SED_wavelen, model_flux_Jy
# %% 

def SED_plot(gaia_id, Teff, mettalicity, log_g, unit):
    SED_wavelen, model_flux_Jy= SED_interpolator(Teff,mettalicity,log_g)
    model_flux = flux_unit_change(model_flux_Jy, unit)

    plt.plot(SED_wavelen, model_flux)
    plt.xlim(0,5)
    plt.xlabel('Wavelength (μm)')
    plt.ylabel('Flux (unit)')
    plt.grid()
    plt.title('SED of stellar model')
    plt.show()

    band_wavelen, flux_values_Jy = get_flux_values(gaia_id)
    flux_values = flux_unit_change(flux_values_Jy, unit)
    flux_val = []
    flux_unc = []
    for i in range(len(band_wavelen)):
        flux_val.append(flux_values[i].value.nominal_value)
        flux_unc.append(flux_values[i].value.std_dev)

    plt.errorbar(band_wavelen, flux_val, yerr=flux_unc, fmt='o',ecolor='orange')
    plt.xlim(0,5)
    plt.xlabel('Wavelength (μm)')
    plt.ylabel('Flux (unit)')
    plt.grid()
    plt.title('SED observed')
    plt.show()

# %% 
gaia_id = 2135550755683407232
Teff = 5785
mettalicity = 0.09
log_g = 4.37

SED_plot(gaia_id, Teff, mettalicity, log_g,'SI')
<<<<<<< HEAD
>>>>>>> update
=======
>>>>>>> update
