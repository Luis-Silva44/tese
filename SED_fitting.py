# %%
from get_flux_values import get_flux_values
from scipy.interpolate import LinearNDInterpolator
import pysynphot as S
import numpy as np
import matplotlib.pyplot as plt

# %% 
first_temp_steps = list(range(3000,13001,250))
second_temp_steps = list(range(13000, 50001, 1000)) 
Teff_grid = first_temp_steps + second_temp_steps
logg_grid = [0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0]
mettalicity_grid = [-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5]

def create_SEDs(Teff_vals, mettalicity_vals, logg_vals):
    SED_data = []

    for teff in Teff_vals:
        for mett in mettalicity_vals:
            for logg in logg_vals:
                #try:
                sed_values = S.Icat('ck04models',teff,mett,logg)
                SED_data.append(((teff,mett,logg),sed_values))
                #except Exception as e:
                #    print(f"Error getting SED values for Teff={teff}, log_g={logg} and mettalicity={mett}")

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
    wavelen = SED_data[0][1].wave
    fluxes = []
    points = []

    for (parameters,SED_values) in SED_data: 
        fluxes.append(SED_values.flux)
        points.append(parameters)
    
    fluxes = np.array(fluxes)
    points = np.array(points)
    
    interpolated_fluxes = []

    for i in range(len(wavelen)):
        flux_interpolator = LinearNDInterpolator(points, fluxes[:,i])
        interpolated_fluxes.append(flux_interpolator(Teff,mettalicity,logg))

    interpolated_fluxes = np.array(interpolated_fluxes)

    wavelen = wavelen * 1e-4
    interpolated_fluxes = interpolated_fluxes * 1e-10 

    return wavelen, interpolated_fluxes

# %% 

def SED_plot(gaia_id, Teff, mettalicity, log_g):
    wavelen, interpolated_fluxes = SED_interpolator(Teff,mettalicity,log_g)
    plt.plot(wavelen, interpolated_fluxes)
    plt.xlim(0,5)
    plt.xlabel('Wavelength (μm)')
    plt.ylabel('Flux (W / cm-2 μm-1)')
    plt.grid()
    plt.title('SED of stellar model')
    plt.show()

    wavelen2, flux_values2, _ = get_flux_values(gaia_id)
    plt.plot(wavelen2, flux_values2, 'o')
    plt.xlim(0,5)
    plt.xlabel('Wavelength (μm)')
    plt.ylabel('Flux (W / cm-2 μm-1)')
    plt.grid()
    plt.title('SED observed')
    plt.show()

# %% 
#gaia_id = '2135550755683407232'
#Teff = 5785
#mettalicity = 0.09
#log_g = 4.37

# SED_plot(gaia_id, Teff, mettalicity, log_g)