# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %% [markdown]
# "SELECT TOP 1000000 gaia_source.source_id,gaia_source.ra,gaia_source.ra_error, \
#                                gaia_source.dec,gaia_source.dec_error,gaia_source.parallax,gaia_source.parallax_error, \
#                                gaia_source.parallax_over_error,gaia_source.phot_g_mean_mag,gaia_source.bp_rp, \
#                                gaia_source.radial_velocity,gaia_source.radial_velocity_error, \
#                                gaia_source.phot_variable_flag,gaia_source.teff_val,gaia_source.a_g_val \
#             FROM gaiadr2.gaia_source \
#             WHERE (gaiadr2.gaia_source.ra>=10.56014098062107 AND gaiadr2.gaia_source.ra<=334.56014098062107 \
#                AND gaiadr2.gaia_source.dec>=-75.587831532659646 AND gaiadr2.gaia_source.dec<=89.587831532659646 \
#                AND gaiadr2.gaia_source.parallax>=11.468624923468156 AND gaiadr2.gaia_source.parallax<=44.068624923468156 \
#                AND gaiadr2.gaia_source.parallax_over_error>=1)"
# Simbad, Ned, Vizier

# %%
import matplotlib.pyplot as plot
import numpy as np
from astropy.io.votable import parse_single_table
import astropy.units as u
import astropy.constants as c
from astropy.table import QTable
import astropy.coordinates as coord
from matplotlib import colors
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D
from PyAstronomy import pyasl
import warnings
warnings.filterwarnings('ignore')
from astroquery.gaia import Gaia
import utilities as util


# %%
#Querying methods

def read_gaia_vot(filename):
    return QTable(parse_single_table(filename).to_table())

def get_table_gaia_query(myquery):
    return QTable(Gaia.launch_job_async(myquery, dump_to_file = True).get_results())


# %%
#Astronomy methods

def luminosityFromMag(magnitude):
    return ((c.L_bol0).value * 10 ** (-0.4 * (magnitude).value))*u.erg

def distance(parallax):
    parallax_arcsec = parallax.to('arcsec').value
    return 1/parallax_arcsec * u.pc

def Magnitude_absolute_from_apparent(mag,dist):
    M = mag - 5*np.log10(dist.to('pc')/u.pc)*u.mag + 5.*u.mag
    return M

def angle_of_star(x,y,z):
    l=np.arctan2(y,x)
    b=np.arctan(z/np.sqrt(x**2+y**2))
    return l,b


# %%
a="trappist-1e"


# %%
ra, dec, pmra, pmdec = util.exoplanetProperMotion(a)
raminor = np.maximum(ra-350, 1)
ramajor = np.minimum(ra+350, 359)
decminor = np.maximum(dec-170, -89)
decmajor = np.minimum(dec+170, 89)


# %%
myquery1 = "SELECT TOP 500 *            FROM gaiadr2.gaia_source             WHERE CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),CIRCLE('ICRS',COORD1(EPOCH_PROP_POS("+ str(ra) +","+ str(dec) +",0,"+ str(pmra) +","+ str(pmdec) +",0,2000,2015.5)),COORD2(EPOCH_PROP_POS("+ str(ra) +","+ str(dec) +",0,"+ str(pmra) +","+ str(pmdec) +",0,2000,2015.5)),0.001388888888888889))=1"

b = get_table_gaia_query(myquery1)

for row in b:
    plx=row['parallax']
print(plx)
plxmajor=(np.minimum((plx.value)+20, 89))
plxminor=(np.maximum((plx.value)-70, 1))


# %%
#Getting the stars between 1.4 and 2 parsecs away and ra between 281 and 310, and dec between 35 and 50

myquery2 = "SELECT TOP 1000000 gaia_source.source_id,gaia_source.ra,gaia_source.ra_error,                                gaia_source.dec,gaia_source.dec_error,gaia_source.parallax,gaia_source.parallax_error,                                gaia_source.parallax_over_error,gaia_source.phot_g_mean_mag,gaia_source.bp_rp,                                gaia_source.radial_velocity,gaia_source.radial_velocity_error,                                gaia_source.phot_variable_flag,gaia_source.teff_val,gaia_source.a_g_val, gaia_source.astrometric_pseudo_colour             FROM gaiadr2.gaia_source             WHERE (gaiadr2.gaia_source.ra>="+ str(raminor) + " AND gaiadr2.gaia_source.ra<="+ str(ramajor) +"               AND gaiadr2.gaia_source.dec>="+ str(decminor) +" AND gaiadr2.gaia_source.dec<="+ str(decmajor) +"                AND gaiadr2.gaia_source.parallax>="+ str(plxminor) +" AND gaiadr2.gaia_source.parallax<="+ str(plxmajor) +"                AND gaiadr2.gaia_source.parallax_over_error>=1)"

p = get_table_gaia_query(myquery2)

p


# %%
largedistance = distance(p['parallax'])
K2Sky=coord.SkyCoord(ra=p['ra'],dec=p['dec'],distance=largedistance)
K2GC=K2Sky.galactocentric
K2Sky


# %%
#Finding Sky Coordinates of Kepler 452
ra= ra*u.degree
dec= dec*u.degree
distance=distance(plx)
K218b=coord.SkyCoord(ra,dec,distance)
K218b


# %%
#Converting Kepler 4452 Sky Coordinates into Galactocentric Coordinates
K218bG=K218b.galactocentric
K218bG


# %%
#Finding the distances between the x,y,z of the stars in our query and Kepler 452's x,y,z
x=(K2GC.x-K218bG.x)
y=(K2GC.y-K218bG.y)
z=(K2GC.z-K218bG.z)

print(x)
print(y)
print(z)


# %%
#Using our method angle_of_star to solve for L,B in Galactic
L,B=angle_of_star(x,y,z)


# %%
#Finding the apparent magnitudes of the stars from the neighboring stars
p['Mg'] = Magnitude_absolute_from_apparent(p['phot_g_mean_mag'], largedistance)
p['otherDistance'] = np.sqrt(x.value**2+y.value**2+z.value**2)*u.pc
p['m'] = (p['Mg'].value - 5 + 5*np.log10(p['otherDistance'].value))*u.mag
p['wavelength']=((1/p['astrometric_pseudo_colour'])*1000)*u.nm/u.um
p['wavelength']


# %%
fig = plot.figure(figsize=(12.8, 9.6))
ax = fig.add_subplot(1,1,1)


plot.scatter((L.to('degree')).value, (B.to('degree').value), s = 1, c = (p['wavelength']), alpha = 0.1, cmap = "nipy_spectral")
ax.set_facecolor("Black")
plot.clim(400, 800)

plot.scatter((L.to('degree')).value, (B.to('degree').value), s = 1, c = (p['m']), alpha = 0.1, cmap = "gray")
ax.set_facecolor("Black")

plot.xlabel("Right Ascension (Degrees)")
plot.ylabel("Declination (Degrees)")

plot.xlim(180,-180)
plot.ylim(-90,90)
plot.show()


# %%
len(p)


# %%



