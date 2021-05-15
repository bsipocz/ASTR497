#!/usr/bin/env python
# coding: utf-8

# # <section class="challenge panel panel-success">
# <div class="panel-heading">
# <h2><span class="fa fa-pencil"></span> Challenge</h2>
# </div>
# 
# 
# <div class="panel-body">
# 
# <ol>
# <li>Examine the headers of the first and second HDU in the point source catalog and try adding new keywords to them or changing them</li>
# <li>Make a histogram of the G-band magnitude of the sources in the catalog - can you figure out what the upper limit for Gmag was when the table was selected?</li>
# <li>Make a plot of the position of the sources on the sky in the point source catalog, in Galactic coordinates l and b (you will need to use what we learned about SkyCoord in the previous tutorial)</li>
# <li>Try and produce a new FITS file that contains the image in the primary HDU and the above table in the second HDU</li>
# </ol>
# 
# </div>
# 
# </section>
# 

# In[1]:


get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
plt.rc('image', origin='lower')
plt.rc('figure', figsize=(10, 6))


# In[2]:


import astropy.units as u
import numpy as np
from astropy.io import fits
from astropy.coordinates import *
from astropy import visualization
from astropy.table import Table


# In[3]:


hdulist1 = fits.open('data/LMCDensFits1k.fits')
hdulist2 = fits.open('data/gaia_lmc_psc.fits')


# In[4]:


hdulist1[0].header


# In[5]:


hdulist2[1].header


# In[6]:


hdulist2[1].header['BITPIX'] = 'YOLO'
hdulist2[1].header['BLARG'] = 'does this work?'
hdulist2[1].header


# In[7]:


hdu2 = Table.read('data/gaia_lmc_psc.fits')
hdu2


# In[8]:


_ = plt.hist(hdu2['phot_g_mean_mag'], bins=100)


# In[9]:


coord1 = SkyCoord(hdu2['ra']*u.deg, hdu2['dec']*u.deg)
coord2 = coord1.transform_to('galactic')
coord2
plt.scatter(coord2.l,coord2.b)


# In[10]:


hdu = fits.PrimaryHDU()
hdu


# In[11]:


hdu.table = hdu2
dataz = fits.getdata('data/gaia_lmc_psc.fits')
hdu2.data = dataz


# In[12]:


hdu.writeto('newFits1.fits', overwrite=True)


# # <section class="challenge panel panel-success">
# <div class="panel-heading">
# <h2><span class="fa fa-pencil"></span> Challenge</h2>
# </div>
# 
# 
# <div class="panel-body">
# 
# <ol>
# <li>Make a table that contains three columns: <code>spectral type</code>, <code>temperature</code>, and <code>radius</code>, and incude 5 rows with fake data (or real data if you like, for example from <a href="http://www.atlasoftheuniverse.com/startype.html">here</a>). Try including units on the columns that can have them.</li>
# <li>Find the mean temperature and the maximum radius</li>
# <li>Try and find out how to add and remove rows</li>
# <li>Add a new column which gives the luminosity (using $L=4\pi R^2 \sigma T^4$)</li>
# </ol>
# 
# </div>
# 
# </section>
# 

# In[13]:


import astropy.constants as const


# In[14]:


table1 = Table()
table1['Spectral Type'] = ['A','B','C','D','E']
table1['Temperature'] = [100,2213,321,4131,5241]*u.K
table1['Radius'] = [5, 21, 2, 34, 12]*const.R_sun
table1


# In[15]:


np.mean(table1['Temperature']), np.max(table1['Radius'])


# In[16]:


table1.add_row(['F', 798*u.K, 6])
table1


# In[17]:


table1['Luminosity']=4*np.pi*const.sigma_sb*[table1['Radius']**2+table1['Temperature']**4]
table1


# # <section class="challenge panel panel-success">
# <div class="panel-heading">
# <h2><span class="fa fa-pencil"></span> Challenge</h2>
# </div>
# 
# 
# <div class="panel-body">
# 
# <p>Starting from the <code>obs</code> table:</p>
# <ol>
# <li>Make a new table that shows every other row, starting with the second row? (that is, the second, fourth, sixth, etc. rows).</li>
# <li>Make a new table the only contains rows where <code>name</code> is <code>M31</code></li>
# </ol>
# 
# </div>
# 
# </section>
# 

# In[18]:


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


# In[19]:


obs2 = obs[2::2]
obs2


# In[20]:


mask1 = obs['name']=='M31'
obs3 = obs[mask1]
obs3


# # <section class="challenge panel panel-success">
# <div class="panel-heading">
# <h2><span class="fa fa-pencil"></span> Challenge</h2>
# </div>
# 
# 
# <div class="panel-body">
# 
# <p>Using the <code>t6</code> table above:</p>
# <ol>
# <li>
# <p>Make a plot that shows <code>j_m</code>-<code>h_m</code> on the x-axis, and <code>h_m</code>-<code>k_m</code> on the y-axis</p>
# </li>
# <li>
# <p>Make a new table that contains the subset of rows where the <code>j_snr</code>, <code>h_snr</code>, and <code>k_snr</code> columns, which give the signal-to-noise-ratio in the J, H, and K band, are greater than 10, and try and show these points in red in the plot you just made.</p>
# </li>
# <li>
# <p>Make a new table (based on the full table) that contains only the RA, Dec, and the <code>j_m</code>, <code>h_m</code> and <code>k_m</code> columns, then try and write out this catalog into a format that you can read into another software package. For example, try and write out the catalog into CSV format, then read it into a spreadsheet software package (e.g. Excel, Google Docs, Numbers, OpenOffice). You may run into an issue at this point - if so, take a look at https://github.com/astropy/astropy/issues/7357 to see how to fix it.</p>
# </li>
# </ol>
# 
# </div>
# 
# </section>
# 

# In[21]:


from astropy.time import Time
from astropy.coordinates import SkyCoord


# In[22]:


#t6 = Table()


# In[23]:


#t6['time'] = Time([50000, 51000, 52000], format='mjd')


# In[24]:


#t6['coord'] = SkyCoord([1, 2, 3] * u.deg, [4, 5, 6] * u.deg)


# In[25]:


#t6['flux'] = [1, 5, 4] * u.mJy


# In[26]:


t6 = Table.read('data/2mass.tbl', format='ascii.ipac')


# In[27]:


subset1 = t6['j_m']-t6['h_m']
subset2 = t6['h_m']-t6['k_m']
plt.scatter(subset1,subset2)


# In[28]:


greaterz = 10
masker = (t6['j_snr'] > greaterz) & (t6['h_snr'] > greaterz) & (t6['k_snr'] > greaterz)
sigToNoise = t6[masker]
plt.scatter(subset1,subset2)
plt.scatter(sigToNoise['j_m']-sigToNoise['h_m'],sigToNoise['h_m']-sigToNoise['k_m'], c='red')


# In[29]:


newTable = Table()
newTable['ra'], newTable['dec'], newTable['j_m'], newTable['h_m'], newTable['k_m'] = t6['ra'], t6['dec'], t6['j_m'], t6['h_m'], t6['k_m']
newTable


# In[30]:


newTable.write('testnewTable.csv', format='csv', overwrite=True)

