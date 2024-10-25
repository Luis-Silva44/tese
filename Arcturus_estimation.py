# %%
from SED_fitting import * 
from get_flux_values import * 
from astroquery.simbad import Simbad
from astroquery.vizier import Vizier
import astropy.units as u
from astropy.coordinates import Angle

# %%  
Teff_arcturus = 4286
logg_Arcturus = 1.66
mettalicity_Arcturus = -0.52

Arcturus_wavelen, Arcturus_interp_fluxes = SED_interpolator(Teff_arcturus, mettalicity_Arcturus, logg_Arcturus)

plt.plot(Arcturus_wavelen, Arcturus_interp_fluxes)
plt.xlim(0,5.0)

# %% 
#print(mag_to_flux(4.915, 'W1'))
#print(mag_to_flux(4.783, 'W2'))
print(mag_to_flux(-2.52, 'J'))
print(mag_to_flux(-2.810, 'H'))
print(mag_to_flux(-2.911, 'K'))

arcturus_fluxes = {'G':8.45e-15, 'GBP':1.07e-14, 'GRP':6.53e-15, 
                   'W1':8.846e-17, 'W2':2.949e-17,
                   'J':2.48e-15, 'H': 1.06e-15, 'K': 4.63e-16}

filter_bands = {'G': 0.673,'GBP':0.532, 'GRP':0.797, 
                        'W1':3.4, 'W2':4.6, 
                        'J':1.25, 'H':1.65, 'K':2.15,}

filter_wavelen = np.array([d for d in filter_bands.values()])
arcturus_flux = np.array([f for f in arcturus_fluxes.values()])

print(filter_wavelen)
print(arcturus_flux)

plt.plot(filter_wavelen, arcturus_flux, 'o')