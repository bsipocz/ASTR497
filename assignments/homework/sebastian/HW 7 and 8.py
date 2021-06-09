#Notebook 7
hdu_gaia = fits.open('data/gaia_lmc_psc.fits')
hdu_gaia.info()

hdu_gaia[1].header

hdu_gaia_ra = hdu_gaia[1].data['ra']
hdu_gaia_dec = hdu_gaia[1].data['dec']
hdu_gaia_ra, hdu_gaia_dec

wcs_gaia = WCS(hdu_gaia[0].header)
wcs_pix = wcs_gaia.wcs_world2pix(hdu_gaia_ra, hdu_gaia_dec, 0)
wcs_pix

import astropy.units as u
catalog = SkyCoord(ra=hdu_gaia_ra * u.deg, dec=hdu_gaia_dec * u.deg)
catalog

gaia_px, gaia_py = wcs.world_to_pixel(catalog)

ax = plt.subplot()
ax.imshow(hdulist[0].data)
ax.plot(gaia_px, gaia_py, '.', color = 'gold')

#------------------------------

#notebook 8

from scipy.ndimage import gaussian_filter
from scipy import misc
import astropy.units as u
from astropy.coordinates import SkyCoord

hdulist_gaia = fits.open('data/gaia_lmc_psc.fits')

hdulist_gauss = gaussian_filter(hdulist[0], sigma=3)
sqrt_norm = simple_norm(hdulist[0].data, stretch='sqrt', percent = 99.5)

plt.imshow(hdulist_iras[0].data, norm = sqrt_norm)
#check to see what the base image looks like first

sqrt_ira = simple_norm(hdulist_iras[0].data, stretch='sqrt', percent = 99.5)
log_ira = simple_norm(hdulist_iras[0].data, stretch='log', percent = 99.5)

ax = plt.subplot(projection = wcs_iras)

plt.imshow(hdulist_iras[0].data, norm = sqrt_norm)
ax.contour(gaussian_filter(hdulist_iras[0].data, sigma=3), 
           transform=ax.get_transform(wcs_iras), colors = '#EA0B07',
           levels = [50, 100, 250, 500])

hdu_g_ra = hdulist_gaia[1].data['ra']
hdu_g_dec = hdulist_gaia[1].data['dec']
hdu_g_ra, hdu_g_dec

wcs_gaia = WCS(hdulist_gaia[0].header)

catalog = SkyCoord(ra=hdu_g_ra * u.deg, dec=hdu_g_dec * u.deg)
catalog

gaia_px, gaia_py = wcs.world_to_pixel(catalog)

ax = plt.subplot()
ax.imshow(hdulist[0].data, norm = sqrt_norm)
ax.plot(gaia_px, gaia_py, '.', color = 'tomato')