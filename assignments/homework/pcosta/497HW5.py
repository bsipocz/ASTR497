from astropy.modeling import fitting
from astropy.modeling import models
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from photutils import CircularAperture
from photutils import aperture_photometry
from photutils import DAOStarFinder
from astropy.stats import sigma_clipped_stats
from photutils.datasets import load_star_image
from astropy.stats import SigmaClip
from photutils import Background2D, MedianBackground
from astropy.table import Table

#Modeling Notebook
#1
x = np.linspace(-5., 5., 200)
y = 3 * np.exp(-0.5 * (x - 1.3)**2 / 0.8**2)
y += np.random.normal(0., 0.2, x.shape)
yerr = np.random.uniform(0.15, 0.25, x.shape)
fitter = fitting.LevMarLSQFitter()

g_init = models.Gaussian1D()
l_init = models.Lorentz1D()
t_init = models.Trapezoid1D()

g_fit_witherr = fitter(g_init, x, y, weights=1/yerr)
l_fit_witherr = fitter(l_init, x, y, weights=1/yerr)
t_fit_witherr = fitter(t_init, x, y, weights=1/yerr)

plt.errorbar(x, y, yerr=yerr, fmt='.', lw=1)
xfine = np.linspace(-5, 5, 1000)
plt.plot(xfine, g_fit_witherr(xfine), lw=2, color='orange')
plt.plot(xfine, l_fit_witherr(xfine), lw=2, color='green')
plt.plot(xfine, t_fit_witherr(xfine), lw=2, color='red')
plt.show()

#2
hdu_gaia = fits.open('../../../astropy_notebooks/data/LMCDensFits1k.fits')[0]
x2, y2 = np.mgrid[:750, :1000]
z = hdu_gaia.data
gg_init = models.Gaussian2D()
gg_fit = fitter(gg_init, x2, y2, z)
plt.imshow(gg_fit(x2, y2))
plt.show()

#3
y_mod = y +0.05 * x**2 + 0.2 * x + 2
combined_init = models.Gaussian1D() + models.Polynomial1D(degree=2)
combined_fit = fitter(combined_init, x, y_mod, weights=1/yerr)
plt.errorbar(x, y_mod, yerr=yerr, fmt='.', lw=1)
plt.plot(xfine, combined_fit(xfine), lw=2)
plt.plot(xfine, combined_fit[0](xfine), lw=2)
plt.show()

#Photometry Notebook
#1.1
star_image = load_star_image()
sigma_clip = SigmaClip(sigma=3)
bkg_estimator = MedianBackground()
bkg = Background2D(star_image.data, (265, 265), sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
star_image_nobkg = star_image.data - bkg.background
mean, median, std = sigma_clipped_stats(star_image_nobkg, sigma=3)
daofind = DAOStarFinder(fwhm=5, threshold=10*std)
sources = daofind(star_image_nobkg)

r5_apertures = CircularAperture(list(zip(sources['xcentroid'], sources['ycentroid'])), r=5)
r10_apertures = CircularAperture(list(zip(sources['xcentroid'], sources['ycentroid'])), r=10)

apphot_table_5 = aperture_photometry(star_image_nobkg, r5_apertures)
apphot_table_10 = aperture_photometry(star_image_nobkg, r10_apertures)

plt.scatter(apphot_table_5['xcenter'], apphot_table_5['ycenter'], s=apphot_table_5['aperture_sum']/10000, color='blue')
plt.scatter(apphot_table_10['xcenter'], apphot_table_10['ycenter'], s=apphot_table_10['aperture_sum']/10000, color='red')
plt.show()

#1.2
divTable = apphot_table_5['aperture_sum'] / apphot_table_10['aperture_sum']
boolean = (divTable < 2.0) & (divTable > 0.5)
commonTable = sources[boolean]
print(commonTable)



