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