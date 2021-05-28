from astropy import units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.visualization import simple_norm
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt
 
plt.rc('image', origin='lower')
plt.rc('figure', figsize=(10, 6))

# Challenge:
# Find the position of all the stars in the data/gaia_lmc_psc.fits catalog in pixel coordinates, and overplot them on the image (you'll need to use things we have learned in previous tutorials).
hdu_image = fits.open('../../../astropy_notebooks/data/LMCDensFits1k.fits')
wcs_image = WCS(hdu_image[0].header)

hdu_gaia = fits.open('../../../astropy_notebooks/data/gaia_lmc_psc.fits')
stars = SkyCoord(ra=hdu_gaia[1].data['ra'], dec=hdu_gaia[1].data['dec'], unit='deg')
px, py = wcs_image.world_to_pixel(stars)

plt.imshow(hdu_image[0].data);
plt.plot(px, py, 'w*', markersize=2);
plt.xlabel('Galactic Longitude')
plt.ylabel('Galactic Latitude')
plt.grid();
plt.show()


# Challenge:
# 1. Make a figure of the IRAS data used above, with the GAIA source density map shown as a contour (note that you might need to smooth the GAIA source density image - check the scipy.ndimage module for some useful functions!)
hdu_iras = fits.open('../../../astropy_notebooks/data/ISSA_100_LMC.fits')
wcs_iras = WCS(hdu_iras[0].header)
sqrt_norm = simple_norm(hdu_iras[0].data, stretch = 'log', percent=97)
ax = plt.subplot(projection=wcs_iras)
ra = ax.coords[0]
dec = ax.coords[1]
ra.set_axislabel("Right Ascension")
dec.set_axislabel("Declination")
ax.imshow(hdu_iras[0].data, norm = sqrt_norm)
ax.contour(gaussian_filter(hdu_image[0].data, 3), transform=ax.get_transform(wcs_image), colors = 'black');
ax.grid()
ax.set_xlim(-0.5, 499.5)
ax.set_ylim(-0.5, 499.5);
ax.show()

# 2. Add the positions of the GAIA sources from the table used in previous tutorials to the image
ax = plt.subplot(projection=wcs_iras)
ra = ax.coords[0]
dec = ax.coords[1]
ra.set_axislabel("Right Ascension")
dec.set_axislabel("Declination")
ax.imshow(hdu_iras[0].data, norm = sqrt_norm)
ax.contour(gaussian_filter(hdu_image[0].data, 3), transform=ax.get_transform(wcs_image), colors = 'black');
ax.plot(stars.ra, stars.dec, '.w', transform=ax.get_transform('world'))
ax.grid()
ax.set_xlim(-0.5, 499.5)
ax.set_ylim(-0.5, 499.5);
ax.show()

# 3. If you have FITS images available, try this out with your own data!
print('I do not have my FITS images readily available, sorry.')
