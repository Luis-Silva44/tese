#%% 
import astropy.units as u
from astroquery.vizier import Vizier
import numpy as np
 

#STILL NEEDED: CORRECT UNCERTAINTIES; REMOVE WARNING FROM COORDINATES
# %% 

def flux_unit_change(value,unit):
    if unit == 'Jy':
        return value.to(u.Jy)
    
    elif unit == 'cgs':
        cgs_flux_units =  u.erg / u.cm**2 / u.s / u.Hz
        return value.to(cgs_flux_units)
    
    elif unit == 'SI':
        SI_flux_units = u.watt / u.m**2 / u.Hz
        return value.to(SI_flux_units)
    
    else:
        raise ValueError('Unit not recognized by programm')
    



# %% 
def vizier_coords(gaia_id):
    Vizier.ROW_LIMIT = 1  
    gaia_catalog = "I/355/gaiadr3"  # Gaia DR3 catalog

    gaia_values = Vizier.query_constraints(catalog=gaia_catalog, Source=str(gaia_id))
    
    if gaia_values:
        ra, dec = gaia_values[0]['RA_ICRS'][0], gaia_values[0]['DE_ICRS'][0]
        return ra, dec
    else:
        raise ValueError("No Gaia id found.")

#%% 
def wise_values(gaia_id):
    wise_catalog = 'II/311/wise'
    ra, dec = vizier_coords(str(gaia_id))
    wise_values = Vizier.query_region(f"{ra} {dec}", radius=10* u.arcsec , catalog=wise_catalog)
    if wise_values:
        return wise_values
    else:
        raise ValueError('No WISE data found')

def two_mass_values(gaia_id):
    two_mass_catalog = 'II/246/out'
    ra, dec = vizier_coords(str(gaia_id))

    two_mass_values = Vizier.query_region(f"{ra} {dec}", radius=10*u.arcsec, catalog=two_mass_catalog)
    
    if two_mass_values:
        return two_mass_values
    else: 
        raise ValueError('No 2MASS data found')

def gaia_values(gaia_id):
    gaia_catalog = "I/355/gaiadr3"  # Gaia DR3 catalog
    gaia_values = Vizier.query_constraints(catalog=gaia_catalog, Source=str(gaia_id))

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
    return flux_constant * 10**(-mag/0.4)

#%%
def get_flux_values(gaia_id):

    filter_bands = {'G': 0.673,'GBP':0.532, 'GRP':0.797, 
                    'W1':3.4, 'W2':4.6, 
                    'J':1.25, 'H':1.65, 'K':2.15,}

    filter_wavelen = np.array([d for d in filter_bands.values()]) * u.um

    wise_data = wise_values(gaia_id)
    two_mass_data = two_mass_values(gaia_id)
    gaia_data = gaia_values(gaia_id)
<<<<<<< HEAD
    W1_mag, W1_mag_err = float(wise_data[0]['W1mag']), float(wise_data[0]['e_W1mag'])
    W2_mag, W2_mag_err = float(wise_data[0]['W2mag']), float(wise_data[0]['e_W2mag'])
    print(W1_mag, W1_mag_err)
    J_mag, J_mag_err = float(two_mass_data[0]['Jmag']), float(two_mass_data[0]['e_Jmag'])
    H_mag, H_mag_err = float(two_mass_data[0]['Hmag']), float(two_mass_data[0]['e_Hmag'])
    K_mag, K_mag_err = float(two_mass_data[0]['Kmag']), float(two_mass_data[0]['e_Kmag'])
=======

    W1_mag = float(wise_data[0]['W1mag'])
    W2_mag = float(wise_data[0]['W2mag'])
    J_mag = float(two_mass_data[0]['Jmag'])
    H_mag = float(two_mass_data[0]['Hmag'])
    K_mag = float(two_mass_data[0]['Kmag'])
>>>>>>> update

    W1_flux = mag_to_flux(W1_mag,'W1')
    W2_flux= mag_to_flux(W2_mag,'W2')
    J_flux = mag_to_flux(J_mag,'J')
    H_flux = mag_to_flux(H_mag,'H')
    K_flux = mag_to_flux(K_mag,'K')

    G_flux = float(gaia_data[0]['FG'])
    GBP_flux = float(gaia_data[0]['FBP'])
    GRP_flux = float(gaia_data[0]['FRP'])
    
    G_flux = G_flux * 1.346109e-21 * u.watt / u.m**2 / u.nm
    GBP_flux = GBP_flux * 3.009167E-21 * u.watt / u.m**2 / u.nm
    GRP_flux = GRP_flux * 1.638483E-21 * u.watt / u.m**2 / u.nm

    G_flux = G_flux.to(u.watt / u.um / u.cm**2)
    GBP_flux = GBP_flux.to(u.watt / u.um / u.cm**2)
    GRP_flux = GRP_flux.to(u.watt / u.um / u.cm**2)

    flux_values = np.array([G_flux.value,GBP_flux.value,GRP_flux.value, W1_flux,W2_flux,J_flux,H_flux,K_flux]) * u.watt / u.um / u.cm**2
    flux_values_Jy = flux_values.to(u.Jy, equivalencies=u.spectral_density(filter_wavelen))

    return filter_wavelen, flux_values_Jy

# %% 
#wavelen, flux_values= get_flux_values(2135550755683407232)
#flux_cgs = flux_unit_change(flux_values, 'cgs')
#flux_Jy = flux_values
#flux_SI = flux_unit_change(flux_values, 'SI')

<<<<<<< HEAD
plt.plot(wavelen, flux_values,'o')
plt.xlabel('Wavelenght (μm)')
plt.ylabel('Flux (W / cm-2 μm-1)')
plt.title('Flux values of the star for each filter')
#plt.errorbar(wavelen, flux_values, yerr=flux_error)

print(flux_values)
print(flux_error)
=======
#plt.plot(wavelen, flux_SI,'o')
#plt.xlabel('Wavelenght (μm)')
#plt.ylabel('Flux (W / m-2 μm-1)')
#plt.title('Flux values of the star for each filter')


#get_flux_values(2135550755683407232)
>>>>>>> update
