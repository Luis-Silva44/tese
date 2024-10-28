#%%
from get_flux_values import * 
from SED_fitting import * 
import pandas as pd

# %% 
def find_nearest_index(array, value):
        index = (np.abs(array - value)).argmin()
        return index

def get_angular_diameter(gaia_id, Teff, mettalicity, log_g):
    wavelen, observed_flux_values = get_flux_values(gaia_id)
    SED_wavelen, SED_fluxes = SED_interpolator(Teff, mettalicity, log_g)

    nearest_index = []
    for i in range(len(wavelen)):
        nearest_index.append(find_nearest_index(SED_wavelen, wavelen [i]))

    model_flux_values = np.array([SED_fluxes[i].value for i in nearest_index])
    
    angular_diameter = 2 * np.sqrt(observed_flux_values / model_flux_values)
    return wavelen, observed_flux_values, model_flux_values, angular_diameter

# %% 
def create_dataframe(gaia_id, Teff, mettalicity, log_g):
    wavelen, observed_flux_values, model_flux_values, angular_diameter = get_angular_diameter(gaia_id, Teff, mettalicity, log_g)

    flux_units_cgs = u.erg / u.cm**2 / u.s / u.Hz
    flux_units_watt = u.watt / u.m**2 / u.Hz

    
    observed_flux_values_cgs = observed_flux_values.to(flux_units_cgs)
    observed_flux_values_watt = observed_flux_values.to(flux_units_watt)

    print(observed_flux_values_cgs)
    print(observed_flux_values_watt)

    flux_table = pd.DataFrame(observed_flux_values, columns=['Flux Observed'])
    flux_table.insert(1, 'Surface flux (model)', model_flux_values)
    flux_table.insert(2, 'Angular Diameter', angular_diameter.value)
    flux_table.insert(0, 'Filter wavelength', wavelen.value)
    flux_table['Angular Diameter'] = flux_table['Angular Diameter'].apply(lambda x: f"{x:.4e}")
    flux_table['Surface flux (model)'] = flux_table['Surface flux (model)'].apply(lambda x: f"{x:.4e}")
    print(flux_table)

    mean_angular_diameter = np.mean(angular_diameter)
    print(mean_angular_diameter)

# %% 

gaia_id = '1019003226022657920'
Teff = 5543
mettalicity = 0.31
log_g = 4.38

SED_plot(gaia_id, Teff, mettalicity, log_g)
create_dataframe(gaia_id, Teff, mettalicity, log_g)