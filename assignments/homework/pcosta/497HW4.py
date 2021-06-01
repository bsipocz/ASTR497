from astropy.convolution import convolve
import numpy as np
from scipy.ndimage import convolve as scipy_convolve
from reproject import reproject_interp
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from astropy.stats import median_absolute_deviation
from photutils.datasets import load_star_image
from astropy.stats import sigma_clip
from astropy.stats import mad_std

#Convolution Notebook
data = np.array([1, 2, 3, 4, 5])
kernel = np.array([2, 2, 4])
default = convolve(data, kernel)
norm = convolve(data, kernel, normalize_kernel=True)
raw = convolve(data, kernel, normalize_kernel=False)
if np.array_equal(default, norm):
    print("Kernel is auto-normalized. Therefore, setting '''normalize_kernel = False''' changes this behavior")
elif np.array_equal(default, raw):
    print("Kernel is not auto-normalized.")

normKernel = np.array([0.2, 0.2, 0.4])
if np.array_equal(scipy_convolve(data, kernel), scipy_convolve(data, normKernel)):
    print("Scipy's convolve auto-normalizes the kernel.")
else:
    print("Scipy's convolve does not auto-normalize the kernel.")

#Reproject Notebook
#1.
hdu_gaia = fits.open('../../../astropy_notebooks/data/LMCDensFits1k.fits')[0]
hdu_iras = fits.open('../../../astropy_notebooks/data/ISSA_100_LMC.fits')[0]
wcs_gaia = WCS(hdu_gaia.header)
wcs_iras = WCS(hdu_iras.header)
gaia_reprojected, footprint = reproject_interp((hdu_gaia.data, wcs_gaia), wcs_iras, shape_out=hdu_iras.data.shape)
plt.imshow(gaia_reprojected, vmax=100)
plt.show()

#2.
ax = plt.subplot(projection=wcs_iras)
ax.imshow(gaia_reprojected, vmax=100)
plt.show()

#Statistics Notebook
statsData = np.hstack([np.random.normal(3, 1, 2000),
                  np.random.normal(10, 0.5, 2000)])
print("Median absolute deviation:", median_absolute_deviation(statsData))

star_image = load_star_image()
clipped_image = sigma_clip(star_image.data, sigma=2, maxiters=20)
print("Robust standard deviation:", mad_std(clipped_image))

