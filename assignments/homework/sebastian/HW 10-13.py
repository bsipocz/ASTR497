#HW 10----

%matplotlib inline
import matplotlib.pyplot as plt
plt.rc('image', origin='lower')
plt.rc('figure', figsize=(10, 6))

from astropy.convolution import Gaussian2DKernel
from astropy.io import fits
from astropy.convolution import convolve

gaia_map = fits.getdata('../../../astropy_notebooks/data/LMCDensFits1k.fits')

gauss = Gaussian2DKernel(3)

ISSA = fits.getdata('../../../astropy_notebooks/data/ISSA_100_LMC.fits')
plt.imshow(convolve(ISSA, gauss))

plt.imshow(scipy_convolve(ISSA, gauss.array))
plt.imshow(convolve(gaia_map, gauss))

#HW 12----

from photutils.datasets import load_star_image
star_image = load_star_image()

from astropy.stats import sigma_clip
clipped_image = sigma_clip(star_image.data, sigma = 2, maxiters=20)

from astropy.stats import median_absolute_deviation
median_absolute_deviation(clipped_image)
