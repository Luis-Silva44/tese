# %% 
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.irsa import Irsa
from astroquery.simbad import Simbad
import pandas as pd

# %% 

ra = 318.5315620098465
dec = +20.7855341949544

coords = SkyCoord(ra = ra*u.deg, dec = dec * u.deg)

result = Irsa.query_region(coords,catalog='fp_psc',spatial='Cone', radius= 5*u.arcmin)
mass_data = result.to_pandas()

print(mass_data.dec)

# %%
result = Simbad.query_region("ASAS J140748-3945.7")
star_data = result.to_pandas()

print(star_data.DEC)
