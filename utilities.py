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
from astroquery.simbad import Simbad

def exoplanetToCoordinates(a):
    result_table = Simbad.query_object(a)
    for row in result_table:
        ra, dec = pyasl.coordsSexaToDeg(row['RA']+ " " + row['DEC'])
    print(ra, dec)
    
def hello():
    print("Hello")
