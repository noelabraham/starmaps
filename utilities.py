from PyAstronomy import pyasl
from astroquery.simbad import Simbad

def exoplanetToCoordinates(a):
    result_table = Simbad.query_object(a)
    for row in result_table:
        ra, dec = pyasl.coordsSexaToDeg(row['RA']+ " " + row['DEC'])
    return ra, dec

def exoplanetProperMotion(a):
    Simbad.add_votable_fields('pm')
    result_table = Simbad.query_object(a)
    for row in result_table:
        ra, dec = pyasl.coordsSexaToDeg(row['RA']+ " " + row['DEC'])
        pmra = row['PMRA']
        pmdec = row['PMDEC']
    return ra, dec, pmra, pmdec


def hello():
    print("Hello")
