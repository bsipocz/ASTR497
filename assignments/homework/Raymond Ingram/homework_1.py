#!/usr/bin/env python
# coding: utf-8

# # HW 1
# Challenge 1: Unit Conversions

# In[1]:

from IPython import get_ipython
get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.units import imperial
import numpy as np


# In[2]:


time = 3 * u.Myr
time


# In[3]:


distance = np.array([1.2, 2.2, 1.7]) * u.kpc
distance


# In[4]:


speed = distance/time
(speed).to(imperial.mile/u.hour)
speed


# In[5]:


pint1 = 1*imperial.pint
pint1


# In[6]:


halfLitre = 0.5*u.l
halfLitre


# In[7]:


if halfLitre < pint1:
    print("A pint is more than a half litre")
else:
    print("A half litre is more than a pint")


# In[8]:


length = 3*u.km
width = 5*u.m
length, width


# In[9]:


area = length*width
area


# In[10]:


area.to(u.m**2)


# In[11]:


area.to(imperial.yd**2)


# Challegne 2: Using equivalencies

# In[12]:


flux = (12*u.mJy).to(u.erg/u.s/u.cm**2/u.Hz)
flux


# ### Challenge 3: Plotting the solar system potential

# In[13]:


from astropy.constants import G, M_sun


# In[14]:


r = np.linspace(1, 50)*u.au
V = -(G*M_sun)/r


# In[15]:


plt.plot(r,V)


# Challenge 4: 

# In[16]:


from astropy.time import Time


# In[17]:


from datetime import datetime
nt = Time.now()
ut = Time(nt, scale='utc', precision=6)
nt, ut


# Challenge 5: SkyCoords

# In[18]:


from astropy.coordinates import SkyCoord, FK5, EarthLocation, AltAz, get_sun


# In[19]:


crabNeb = SkyCoord.from_name('M1').transform_to('galactic');
crabNeb


# In[20]:


tTauri = SkyCoord.from_name('T Tauri');
tTauri


# In[21]:


seper = crabNeb.separation(tTauri).to(u.deg)
seper


# In[22]:


l = np.random.uniform(0, 360, 30) * u.deg
b = np.random.uniform(-1, 1, 30) * u.deg
randPosns = SkyCoord(l, b, frame = 'galactic')
randPosns = randPosns.transform_to(FK5(equinox = 'J2021'))
plt.scatter(randPosns.ra, randPosns.dec, marker='.', color='tab:blue', s=100);
plt.title("Random Positions in the Galactic plane b/w l=0 and l=360 and b=+-1")
plt.xlabel("l, deg");
plt.ylabel("b, deg");


# In[23]:


seattle = EarthLocation(lat = 47.65361111*u.deg, lon = -122.31138889*u.deg, height = 0*u.m)
timeZoneOffset = -7*u.hour  #Sets the time in the graph to PDT (change 7 to 8 during the rest of the year)
midnight = Time('2021-4-29 00:00:00') - timeZoneOffset
timeStep = np.linspace(-12, 12, 1000)*u.hour
times = midnight + timeStep
sun_altaz = get_sun(times).transform_to(AltAz(obstime=times, location=seattle))
plt.scatter(timeStep, sun_altaz.alt, label='Sun', s=8)
plt.fill_between(timeStep, 0, 90, sun_altaz.alt < -0*u.deg, color='0.5', zorder=0, label='Dawn/Dusk')
plt.fill_between(timeStep, 0, 90, sun_altaz.alt < -18*u.deg, color='k', zorder=0, label='Night')
plt.colorbar().set_label('Azimuth [deg]')
plt.legend(loc='upper left')
plt.xlim(-12, 12)
plt.xticks((np.arange(13)*2-12))
plt.title('Sun in Sky, UTC')
plt.ylim(0, 90)
plt.xlabel('Hours from PDT Midnight')
plt.ylabel('Altitude [deg]')


# In[ ]:




