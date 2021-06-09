# imports
import numpy as np
from astropy import units as u
from astropy.modeling import fitting, models
from astropy.io import fits
from astropy.coordinates import SkyCoord
from photutils.datasets import load_star_image
from astropy.stats import sigma_clip, SigmaClip, sigma_clipped_stats
from photutils import CircularAperture, aperture_photometry, DAOStarFinder, Background2D, MedianBackground

import matplotlib.pyplot as plt
plt.rc('image', origin='lower')
plt.rc('figure', figsize=(12, 8))

#solutions
# Challenge: Try fitting a Lorentzian and a Trapezoidal model to the same data. You can use tab-completion to find these models, or search in the documentation page. Overplot these models on the data along with the Gaussian fit.
x1 = np.linspace(-5., 5., 200)
y1 = 3 * np.exp(-0.5 * (x1 - 1.3)**2 / 0.8**2)
y1 += np.random.normal(0., 0.2, x1.shape)
yerr = np.random.uniform(0.15, 0.25, x1.shape)

# fitter
fitter = fitting.LevMarLSQFitter()
# models
g_init = models.Gaussian1D()
l_init = models.Lorentz1D()
t_init = models.Trapezoid1D()

# fit
g_fit = fitter(g_init, x1, y1)
l_fit = fitter(l_init, x1, y1)
t_fit = fitter(t_init, x1, y1)

# errors
g_fit_witherr = fitter(g_init, x1, y1, weights=1/yerr)
l_fit_witherr = fitter(l_init, x1, y1, weights=1/yerr)
t_fit_witherr = fitter(t_init, x1, y1, weights=1/yerr)

# plotting
plt.errorbar(x1, y1, yerr=yerr, fmt='.')
xfine = np.linspace(-5, 5, 1000)
plt.plot(xfine, g_init(xfine), lw=2, label = 'g_init', c='g')
plt.plot(xfine, g_fit(xfine), lw=5, label = 'g_fit', c='b')
plt.plot(xfine, g_fit_witherr(xfine), lw=2, label = 'g_fit w/ err', c='magenta')
plt.plot(xfine, l_init(xfine), lw=2, label = 'l_init', c='r')
plt.plot(xfine, l_fit(xfine), lw=2, label = 'l_fit', c='cyan')
plt.plot(xfine, l_fit_witherr(xfine), lw=2, label = 'l_fit w/ err', c='brown');
plt.plot(xfine, t_init(xfine), lw=2, label = 't_init')
plt.plot(xfine, t_fit(xfine), lw=5, label = 't_fit')
plt.plot(xfine, t_fit_witherr(xfine), lw=2, label = 't_fit w/ err')
plt.title("Fitting Lorentzian and Trapezoidal models", fontsize=20);
plt.legend();
plt.show()


# Challenge
# Read in the GAIA source density FITS file we used in previous tutorials, and try fitting a 2D Gaussian to it.
hdu_gaia = fits.open('astropy_notebooks/data/LMCDensFits1k.fits')[0]
plt.imshow(hdu_gaia.data);
plt.title("GAIA Data", fontsize=20)
plt.show()
p_init = models.Gaussian2D()
y, x = np.mgrid[:750, :1000]
p_fit = fitter(p_init, x, y, z=hdu_gaia.data)
plt.imshow(p_fit(x, y));
plt.title("2D Gaussian fit", fontsize=20)
plt.show()
plt.imshow(hdu_gaia.data - p_fit(x, y));
plt.title("GAIA Data -  2D Gaussian fit", fontsize=20);
plt.show()


# Challenge
# 1. Modify the plot of the combined fit to show the individual model components for the best-fit parameters.
y_mod = y1 + 0.05 * x1**2 +0.2 * x1 + 2
_ = plt.errorbar(x1, y_mod, yerr=yerr, fmt='.', c='blue')

# fitting
combined_init = models.Gaussian1D() + models.Polynomial1D(degree=2)
combined_fit = fitter(combined_init, x1, y_mod, weights=1/yerr)

plt.plot(xfine, combined_fit(xfine), lw=4, c='green', label='combined')
plt.plot(xfine, combined_fit[0](xfine), lw=4, c='red', label='gaussian')
plt.plot(xfine, combined_fit[1](xfine), lw=4, c='magenta', label='polynomial')
plt.title("Combined fits separated", fontsize=20);
plt.show()

# 2. Continue the previous challenge to fit the LMC source density map by fitting two Gaussians to it. Be aware 
# that especially with compound models, initial values matter! [advanced]
print("I'm not sure how to set up the initial values, but the solution would be a combination of the two previous questions.")

# Challenge
# 1. Carry out aperture photometry using a radius of 5 and 10 pixels, and compare the results in a scatter plot
star_image = load_star_image()
sigma_clip = SigmaClip(sigma=3)
bkg_estimator = MedianBackground()
bkg = Background2D(star_image.data, (265, 265), sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
star_image_nobkg = star_image.data - bkg.background
mean, median, std = sigma_clipped_stats(star_image_nobkg, sigma=3)
daofind = DAOStarFinder(fwhm=5, threshold=10*std)
sources = daofind(star_image_nobkg)
ap_5 = CircularAperture(list(zip(sources['xcentroid'], sources['ycentroid'])), r=5)
ap_10 = CircularAperture(list(zip(sources['xcentroid'], sources['ycentroid'])), r=10)
phot_5 = aperture_photometry(star_image_nobkg, ap_5)
phot_10 = aperture_photometry(star_image_nobkg, ap_10)
fig, ax = plt.subplots(2, figsize=(10, 15))
ax[0].scatter(phot_5['xcenter'], phot_5['ycenter'], s=phot_5['aperture_sum']/10000, c='blue')
ax[0].set_title("Photometry with aperture radius 5", fontsize=20);
ax[1].scatter(phot_10['xcenter'], phot_10['ycenter'], s=phot_10['aperture_sum']/10000, c='red');
ax[1].set_title("Photometry with aperture radius 10", fontsize=20)
fig.show()

# 2. Using this, construct a table containing just the sources where the 5 and 10 pixel fluxes agree to within a factor of 2
flux_5 = phot_5['aperture_sum']
flux_10 = phot_10['aperture_sum']
tendiv5= flux_10/flux_5 < 2
fivediv10 = flux_5/flux_10 < 2
flux_factor2 = sources[tendiv5 & fivediv10]
print(flux_factor2)