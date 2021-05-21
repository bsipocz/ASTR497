from astropy import units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.visualization import simple_norm
import matplotlib.pyplot as plt
plt.rc('image', origin='lower')
plt.rc('figure', figsize=(10, 6))

# Challenge:
# Find the position of all the stars in the data/gaia_lmc_psc.fits catalog in pixel coordinates, and overplot them on the image (you'll need to use things we have learned in previous tutorials).
hdu_image = fits.open('astropy_notebooks/data/LMCDensFits1k.fits')
wcs_image = WCS(hdu_image[0].header)
hdu_gaia = fits.open('astropy_notebooks/data/gaia_lmc_psc.fits')
wcs_gaia = WCS(hdu_gaia[1].header)
ra = hdu_gaia[1].data['ra']
dec = hdu_gaia[1].data['dec']
stars = SkyCoord(ra=ra * u.deg, dec=dec * u.deg)
overplot = wcs_image.world_to_pixel(stars)
ax = plt.subplot(projection=wcs_image)
lon = ax.coords[0]
lat = ax.coords[1]
lon.set_axislabel("Galactic Longitude")
lat.set_axislabel("Galactic Latitude")
ax.imshow(hdu_image[0].data);
ax.plot(overplot[0], overplot[1], 'w*', markersize=2);
ax.grid()
plt.show()


# Challenge:
# 1. Make a figure of the IRAS data used above, with the GAIA source density map shown as a contour (note that you might need to smooth the GAIA source density image - check the scipy.ndimage module for some useful functions!)
hdu_iras = fits.open('astropy_notebooks/data/ISSA_100_LMC.fits')
wcs_iras = WCS(hdu_iras[0].header)
sqrt_norm = simple_norm(hdu_iras[0].data, stretch = 'log', percent=97)
ax = plt.subplot(projection=wcs_iras)
ra = ax.coords[0]
dec = ax.coords[1]
ra.set_axislabel("Right Ascension")
dec.set_axislabel("Declination")
ax.imshow(hdu_iras[0].data, norm = sqrt_norm);
ax.contour(hdu_image[0].data, transform=ax.get_transform(wcs_image), colors = 'black', levels=[50, 100, 250, 500]);
ax.grid()
ax.set_xlim(-0.5, 501);
ax.set_ylim(-0.5, 502);
plt.show()

# 2. Add the positions of the GAIA sources from the table used in previous tutorials to the image
ax = plt.subplot(projection=wcs_iras)
ra = ax.coords[0]
dec = ax.coords[1]
ra.set_axislabel("Right Ascension")
dec.set_axislabel("Declination")
ax.imshow(hdu_iras[0].data, norm = sqrt_norm);
ax.contour(hdu_image[0].data, transform=ax.get_transform(wcs_image), colors = 'black', levels=[50, 100, 250, 500]);
ax.set_xlim(-0.5, 501);
ax.set_ylim(-0.5, 502);
ax.plot(overplot[0], overplot[1], 'wo', markersize=2);
ax.grid()
plt.show()

# 3. If you have FITS images available, try this out with your own data!
print('I do not have my FITS images readily available, sorry.')