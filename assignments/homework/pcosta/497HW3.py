from astropy.wcs import WCS
from astropy.io import fits
import matplotlib.pyplot as plt

#WCS Notebook

imagehdulist = fits.open('../../../astropy_notebooks/data/LMCDensFits1k.fits')
plt.imshow(imagehdulist[0].data)
gaiaHDUlist = fits.open('../../../astropy_notebooks/data/gaia_lmc_psc.fits')
gaiaData = gaiaHDUlist[1].data
wcs = WCS('../../../astropy_notebooks/data/LMCDensFits1k.fits')
px, py = wcs.wcs_world2pix(gaiaData['ra'], gaiaData['dec'], 0)
plt.plot(px, py, '.')
plt.show()

#WCSAxes Notebook
#1
LMChdulist = fits.open('../../../astropy_notebooks/data/LMCDensFits1k.fits')
IRAShdulist = fits.open('../../../astropy_notebooks/data/ISSA_100_LMC.fits')
wcs = WCS(IRAShdulist[0].header)
ax = plt.subplot(projection=wcs)
ax.imshow(IRAShdulist[0].data)
ax.grid()
ax.contour(LMChdulist[0].data, transform=ax.get_transform(wcs),
           colors='white', levels=[50,100,250,500])


#2.
plt.plot(gaiaData['ra'], gaiaData['dec'], '.')
plt.show()