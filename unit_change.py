# %% 

import astropy.units as u

# %% 

def flux_unit_change(value,unit):
    if unit == 'cgs':
        cgs_flux_units =  u.erg / u.cm**2 / u.s / u.Hz
        return value.to(cgs_flux_units)
    if unit == 'Jansky':
        return value.to(Jansky)
    if unit == 'SI':
        SI_flux_units = flux_units_watt = u.watt / u.m**2 / u.Hz
        return value.to(SI_flux_units)