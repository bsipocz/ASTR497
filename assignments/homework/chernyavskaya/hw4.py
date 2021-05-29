# imports
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.convolution import Gaussian2DKernel, convolve 
from scipy.ndimage import convolve as scipy_convolve
from reproject import reproject_interp as r_i, reproject_exact as r_e
from photutils.datasets import load_star_image
from astropy.stats import sigma_clip, sigma_clipped_stats, mad_std, median_absolute_deviation, jackknife_stats


import matplotlib.pyplot as plt
plt.rc('image', origin='lower')
plt.rc('figure', figsize=(10, 6))
# solutions

# Challenge
# Using a simple 1D dataset as done above, can you determine whether the kernel is automatically normalized by 
# default? How can you change this behavior? And how does this compare to SciPy's convolve function?
data = [1, 2, np.nan, 4, 5]
kernel = [0.5, 1, 0.5]
apy_conv = convolve(data, kernel)
apy_conv1 = convolve(data, kernel, normalize_kernel=False)
apy_conv2 = convolve(data, kernel)
scipy_conv = scipy_convolve(data, kernel)
print(f"Astropy convolve, default is normalized: {apy_conv} \nAstropy convolve, normalize_kernel==False: {apy_conv1} \nScipy convolve lacks normalization: {scipy_conv}")


# Challenge
# 1. Reproject the GAIA source density map to the WCS of the IRAS map
hdu_gaia = fits.open('../../../astropy_notebooks/data/LMCDensFits1k.fits')[0]
hdu_iras = fits.open('../../../astropy_notebooks/data/ISSA_100_LMC.fits')[0]
wcs_gaia = WCS(hdu_gaia.header)
wcs_iras = WCS(hdu_iras.header)
gaia_reprojected, footprint = r_i(hdu_gaia, hdu_iras.header)
# 2. Make a plot of the resulting image with WCSAxes
#plt.subplot(projection=wcs_gaia)
ax = plt.subplot(projection=wcs_iras)
ax.imshow(gaia_reprojected, vmin = 0, vmax=100);
plt.show()
# 3. If you have FITS images available, try this out with your own data!
print("I don't have FITS files that have WCS info in them.")

# Challenge
# 1. Take a look at the astropy.stats documentation, and in particular the long list of functions at the bottom, in 
# case you see something that could be useful to you! (and feel free to try them if so)
print(f"Something useful to me is the jackknife testing functionality. See below:")
print("estimate, bias, stderr, conf_interval=",jackknife_stats(np.array([1,2,3,4,5,6,7,8,9,0]), np.mean, 0.95))
# 2. If you had to find the median absolute deviation for a dataset, how would you do it? Try and find the robust 
# standard deviation using the median absolute deviation for the sigma clipped array we produced before.
star_image = load_star_image()
clipped_image = sigma_clip(star_image.data, sigma=2, maxiters=20)
mad = median_absolute_deviation(clipped_image)
robust_std_mad = mad_std(clipped_image)
print("\nMedian absolute deviation:", mad, "\nRobust standard deviation using MAD:", robust_std_mad)