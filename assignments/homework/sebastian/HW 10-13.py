#HW 10----

%matplotlib inline
import matplotlib.pyplot as plt
plt.rc('image', origin='lower')
plt.rc('figure', figsize=(10, 6))

from astropy.convolution import Gaussian2DKernel
from astropy.io import fits

gauss = Gaussian2DKernel(3)

ISSA = fits.getdata('data/ISSA_100_LMC.fits')
plt.imshow(convolve(ISSA, gauss))

plt.imshow(scipy_convolve(ISSA, gauss.array))
plt.imshow(convolve(gaia_map, gauss))

#HW 12----

%matplotlib inline
import matplotlib.pyplot as plt
plt.rc('image', origin='lower')
plt.rc('figure', figsize=(10, 6))

from photutils.datasets import load_star_image
star_image = load_star_image()
star_image

from astropy.stats import median_absolute_deviation
median_absolute_deviation(star_image.data)

#HW 13----

from astropy.modeling import fitting

from astropy.io import fits
hdulist = fits.open('data/LMCDensFits1k.fits')

from astropy.wcs import WCS
wcs = WCS(hdulist[0].header)

ax = plt.subplot(projection=wcs)
ax.imshow(hdulist[0].data)

p_gauss = models.Gaussian2D()
p_fit_g = fitter(p_gauss, x2, y2, z2)



plt.errorbar(x, y_mod, yerr=yerr, fmt='.')
plt.plot(xfine, combined_fit(xfine), lw=2)
plt.plot(xfine, combined_fit[0](xfine), lw=2)