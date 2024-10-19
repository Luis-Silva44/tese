#%% 
import astroquery as astro
from astropy import units as u
from astroquery.vizier import Vizier

# %% 

def vizier_coords(gaia_id):
    Vizier.ROW_LIMIT = 1  
    gaia_catalog = "I/355/gaiadr3"  # Gaia DR3 catalog

    gaia_values = Vizier.query_constraints(catalog=gaia_catalog, Source=gaia_id)
    
    if gaia_values:
        ra, dec = gaia_values[0]['RA_ICRS'][0], gaia_values[0]['DE_ICRS'][0]
        return ra, dec
    else:
        print("No Gaia id found.")

#%% 
def wise_values(gaia_id):
    wise_catalog = 'II/311/wise'
    ra, dec = vizier_coords(gaia_id)

    result_wise = Vizier.query_region(f"{ra} {dec}", radius=2* u.arcsec , catalog=wise_catalog)
    
    if result_wise:
        return result_wise
    else:
        print('No WISE data found')

def two_mass_values(gaia_id):
    two_mass_catalog = 'II/246/out'
    ra, dec = vizier_coords(gaia_id)
    two_mass_values = Vizier.query_region(f"{ra} {dec}", radius="5s", catalog=two_mass_catalog)

    if two_mass_values:
        return two_mass_values
    else:
        print("No 2MASS data found.")

def gaia_values(gaia_id):
    gaia_catalog = "I/355/gaiadr3"  # Gaia DR3 catalog
    gaia_values = Vizier.query_constraints(catalog=gaia_catalog, Source=gaia_id)

    if gaia_values:
        return gaia_values
    else:
        print("No Gaia data found.")