import numpy as np
from astropy import units as u
import astropy.constants as const
from astropy.io import fits
from astropy.table import Table, QTable
from astropy.time import Time
from astropy.coordinates import SkyCoord

import matplotlib.pyplot as plt
plt.rc('image', origin='lower')
plt.rc('figure', figsize=(10, 6))

# Challenge
#1. Examine the headers of the first and second HDU in the point source catalog and try adding new keywords to them or changing them
hdu1 = fits.open('../../../astropy_notebooks/data/LMCDensFits1k.fits')
hdu2 = fits.open('../../../astropy_notebooks/data/gaia_lmc_psc.fits')
header1 = hdu1[0].header
header1['telescop'] = 'gaia'
print(header1)
header2 = hdu2[1].header
header2['telescop'] = 'gaia'
print(header2)

#2. Make a histogram of the G-band magnitude of the sources in the catalog - can you figure out what the upper limit for Gmag was when the table was selected?
mags = hdu2[1].data['phot_g_mean_mag']
plt.hist(mags, color= 'green', bins = 25)
plt.xlim(3,12);
plt.ylim(0,300);
plt.xlabel('G band Magnitude', fontsize = 12);
plt.ylabel('Frequency', fontsize = 12);
plt.title('G band Magnitude histogram', fontsize = 20);
plt.show();
print(f'The upper limit of Gmag is {np.max(mags):.3f}.')

#3. Make a plot of the position of the sources on the sky in the point source catalog, in Galactic coordinates l and b (you will need to use what we learned about SkyCoord in the previous tutorial)
ra = hdu2[1].data['ra']
dec = hdu2[1].data['dec']
coords = SkyCoord(ra * u.deg, dec * u.deg)
newcoords = coords.transform_to('galactic')
plt.scatter(newcoords.l, newcoords.b,marker = '*', color = 'blue', edgecolors = 'black', s=100);
plt.xlabel('l', fontsize = 12);
plt.ylabel('b', fontsize = 12);
plt.title('Galactic coordinate plot of sources on sky', fontsize = 20)
plt.show();

#4. Try and produce a new FITS file that contains the image in the primary HDU and the above table in the second HDU
primary_hdu = fits.open('../../../astropy_notebooks/data/LMCDensFits1k.fits')[0]
secondary_hdu = fits.open('../../../astropy_notebooks/data/gaia_lmc_psc.fits')[1]
hdnew = fits.HDUList([primary_hdu, secondary_hdu])
hdnew.writeto('hdnew.fits')


# Challenge
# 1. Make a table that contains three columns: spectral type, temperature, and radius, and incude 5 rows with fake data (or real data if you like, for example from http://www.atlasoftheuniverse.com/startype.html. Try including units on the columns that can have them.
t1 = Table()
t1['spectral type'] = ['O3', 'O5', 'O8', 'B0', 'B3']
t1['temp'] = [53000, 35000, 35000, 30000, 19000] * u.K
t1['radius'] = [15, 12, 8.5, 7.4, 4.8] * const.R_sun
print(t1)

# 2. Find the mean temperature and the maximum radius
print(f"Mean temp is {np.mean(t1['temp']):5.0f}, maximum radius is {np.max(t1['radius']):.2e}.")


# 3.Try and find out how to add and remove rows
t2 = t1
t2.add_row(['B5', 15000*u.K, 3.9*const.R_sun])
print(t2)
t3 = t2
t3.remove_row(5)
print(t3)
print("To add row, use the add_row() method, and to remove a row use remove_row().")

# 4. Add a new column which gives the luminosity (using L=4*pi R^2 *sigma T^4)
t1['luminosity'] = 4 * np.pi * (t1['radius']**2) * const.sigma_sb * (t1['temp']**4)
print(t1)

# Challenge
# Starting from the obs table:
obs = Table(rows=[('M31' , '2012-01-02', 17.0, 17.5),
                  ('M31' , '2012-01-02', 17.1, 17.4),
                  ('M101', '2012-01-02', 15.1, 13.5),
                  ('M82' , '2012-02-14', 16.2, 14.5),
                  ('M31' , '2012-02-14', 16.9, 17.3),
                  ('M82' , '2012-02-14', 15.2, 15.5),
                  ('M101', '2012-02-14', 15.0, 13.6),
                  ('M82' , '2012-03-26', 15.7, 16.5),
                  ('M101', '2012-03-26', 15.1, 13.5),
                  ('M101', '2012-03-26', 14.8, 14.3)],
            names=['name', 'obs_date', 'mag_b', 'mag_v'])
#1. Make a new table that shows every other row, starting with the second row? (that is, the second, fourth, sixth, etc. rows).
evens = obs[1::2]
print(evens)

# 2. Make a new table the only contains rows where name is M31
M31 = obs[obs['name'] == 'M31']
print(M31)

# Challenge
Mass2 = Table.read('astropy_notebooks/data/2mass.tbl', format='ascii.ipac')
# Using the Mass2 table above:
# 1. Make a plot that shows j_m-h_m on the x-axis, and h_m-k_m on the y-axis
x_ax = Mass2['j_m'] - Mass2['h_m']
y_ax = Mass2['h_m'] - Mass2['k_m']
plt.scatter(x_ax, y_ax, marker='*', color='green', s= 75);
plt.xlabel('J - H', fontsize = 12);
plt.ylabel('H - K', fontsize = 12);
plt.title('Colour-colour diagram', fontsize = 20);
plt.show();


# 2. Make a new table that contains the subset of rows where the j_snr, h_snr, and k_snr columns, which give the signal-to-noise-ratio in the J, H, and K band, are greater than 10, and try and show these points in red in the plot you just made.
subset = Mass2[(Mass2['j_snr']>10) & (Mass2['h_snr']>10) & (Mass2['k_snr']>10)]
x_ax = Mass2['j_m'] - Mass2['h_m']
y_ax = Mass2['h_m'] - Mass2['k_m']
x_set = subset['j_m'] - subset['h_m']
y_set = subset['h_m'] - subset['k_m']
plt.scatter(x_ax, y_ax, marker='*', color='k', linewidth = 2, edgecolor ='green', s=150);
plt.scatter(x_set, y_set, marker='*', color = 'k', edgecolor ='red', s=150)
plt.xlabel('J - H', fontsize = 12);
plt.ylabel('H - K', fontsize = 12);
plt.title('Colour-colour diagram', fontsize = 20);
plt.show();


# 3. Make a new table (based on the full table) that contains only the RA, Dec, and the j_m, h_m and k_m columns, then try and write out this catalog into a format that you can read into another software package. For example, try and write out the catalog into CSV format, then read it into a spreadsheet software package (e.g. Excel, Google Docs, Numbers, OpenOffice). You may run into an issue at this point - if so, take a look at https://github.com/astropy/astropy/issues/7357 to see how to fix it.
t2 = Table()
t2 = Mass2['ra', 'dec', 'j_m', 'h_m', 'k_m']

t2.write('t2.csv')
print("Ran into no issues opening in Excel.")