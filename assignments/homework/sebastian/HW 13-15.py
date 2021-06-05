#HW 13----

%matplotlib inline
import matplotlib.pyplot as plt
plt.rc('image', origin='lower')
plt.rc('figure', figsize=(10, 6))

from astropy.modeling import fitting

from astropy.io import fits

hdulist = fits.open('../../../data/LMCDensFits1k.fits')

from astropy.wcs import WCS
wcs = WCS(hdulist[0].header)

ax = plt.subplot(projection=wcs)
ax.imshow(hdulist[0].data)

p_gauss = models.Gaussian2D()
p_fit_g = fitter(p_gauss, x2, y2, z2)



plt.errorbar(x, y_mod, yerr=yerr, fmt='.')
plt.plot(xfine, combined_fit(xfine), lw=2)
plt.plot(xfine, combined_fit[0](xfine), lw=2)

#HW 14----

from photutils import CircularAperture
from photutils import aperture_photometry

aperture5 = CircularAperture(list(zip(sources['xcentroid'], sources['ycentroid'])), r=5)
aperture10 = CircularAperture(list(zip(sources['xcentroid'], sources['ycentroid'])), r=10)

app5_tab = aperture_photometry(star_image_nobkg, aperture5)
app10_tab = aperture_photometry(star_image_nobkg, aperture10)

plt.scatter(app5_tab['xcenter'], app5_tab['ycenter'], s=app5_tab['aperture_sum']/10000)
plt.scatter(app5_tab['xcenter'], app5_tab['ycenter'], s=app5_tab['aperture_sum']/10000)

tab_510 = app5_tab['aperture_sum'] / app10_tab['aperture_sum']

cut_tab = (tab_510 < 2.0) & (tab_510 > 0.5)

filtered = sources[cut_tab]

#HW 15----

from astropy.io import fits
from astropy_healpix import HEALPix
from astropy.coordinates import SkyCoord


hdulist = fits.open('../../../data/HFI_SkyMap_857_2048_R1.10_nominal_ZodiCorrected_lowres.fits')
hdulist.info()

hdulist[1].data['I_STOKES'].mean()

hp.skycoord_to_healpix(SkyCoord.from_name('M42'), hdulist[1].data['I_STOKES'].mean())

hp.cone_search_skycoord(SkyCoord.from_name('M42'), radius=2*u.deg)
