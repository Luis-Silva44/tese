# %% 
from extinction import ccm89, apply
from get_flux_values import gaia_values, get_flux_values
from sympy import solve, symbols
import astropy.units as u
from SED_fitting import SED_interpolator


 #%% 

def intrinsic_color(Teff, mettalicity):
    x = symbols('x')
    sol = solve(8939 - 6395 * x + 2381 * x**2 - Teff + 451 * mettalicity + 154*mettalicity**2)
    if sol[0] < .40 or sol[0] > 1.20:
        sol.pop(0)
    if sol[1] < .40 or sol[1] > 1.20:
        sol.pop(1)   
    return sol

def color_excess(Teff, mettalicity, gaia_id):
    gaia_vals = gaia_values(gaia_id)
    color = float(gaia_vals[0]['BPmag']) - float(gaia_vals[0]['Gmag'])
    int_color =  intrinsic_color(Teff, mettalicity)
    return color - int_color[0]

def flux_extinction(wavelen, flux, Teff, mettalicity, gaia_id):
    col_exc = color_excess(Teff, mettalicity, gaia_id)
    flux_ext = apply(ccm89(wavelen.to(u.angstrom), col_exc*3.1, 3.1), flux)
    return flux_ext

# %% 

Teff = 5431 
mett = -0.45
log_g = 4.33 
gaiaid = 5855730584310531200
wavelen, flux = get_flux_values(gaiaid)
#print(flux)
flux_ext = flux_extinction(wavelen,flux,Teff, mett, 1019003226022657920)
#print(flux_ext)

# %% 
wavelen, flux = SED_interpolator(Teff, mett, log_g)
type(flux)
#flux_ext = apply(ccm89(wavelen.to(u.angstrom), -0.3*3.1, 3.1), flux)