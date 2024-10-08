# %% 
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.irsa import Irsa
from astroquery.simbad import Simbad
import pandas as pd
from astroquery.gaia import Gaia

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
    result = Irsa.query_region(coords, catalog='allwise_p3as_psd', spatial='Cone', radius = .05*u.arcmin)
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
print(wise_designation)

# %% 
two_mass_designation = two_mass_data('20 LMi')
print(two_mass_designation)

# %% 
gaia_designation = GAIA_data("20 LMi")
print(gaia_designation)
