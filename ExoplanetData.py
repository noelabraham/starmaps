from astroquery.nasa_exoplanet_archive import NasaExoplanetArchive

exoplanet_archive_table = NasaExoplanetArchive.get_confirmed_planets_table()

exoplanet_archive_table

exoplanet = NasaExoplanetArchive.query_planet('Trappist-1e')

exoplanet