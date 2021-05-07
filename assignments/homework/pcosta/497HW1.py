import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.units import imperial
from astropy import constants as const
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy.coordinates import EarthLocation, AltAz
from random import *
from astropy.coordinates import FK5

#Units Notebook Homework
#1-Unit Conversions
velocities = np.array([2.3, 3.5, 6.0, 10.5])*u.m/u.s
convertedVelocities = velocities.to(imperial.mi/u.h)
print(convertedVelocities)

#2-Using Equivalencies
specDensity1 = (12 * u.mJy).to(u.erg / u.cm**2 / u.s / u.Hz)
print(specDensity1)

specDensity2 = (12 * u.mJy).to(u.W / u.m**2 / u.Hz)
print(specDensity2)

#3-Plotting the Solar System Potential
def V(r):
    return -1 * const.G * const.M_sun / r
R = ((np.arange(1, 50, 0.5))*u.AU)
plt.plot(R, V(R).to(u.mJ/u.kg))
plt.title("Gravitational Potential")
plt.xlabel("Radius (AU)")
plt.ylabel("Potential (MJ/Kg)")
plt.show()

#Times Notebook Homework
t1 = Time.now()
print(t1)
Start = Time("2021-04-02 12:00:00")
print((t1-Start).to(u.min))

#Coordinates Notebook
#1.
m1 = SkyCoord.from_name('M1')
M1 = m1.transform_to('galactic')
print(M1)

#2.
tauri = SkyCoord.from_name('T Tauri')
seattle = EarthLocation(lat=47.6062*u.deg, lon=-122.3321*u.deg, height=0*u.m)
utc_time = Time.now()
m1_altaz = m1.transform_to(AltAz(obstime=utc_time, location=seattle))
tauri_altaz = tauri.transform_to(AltAz(obstime=utc_time, location=seattle))
separation = (tauri_altaz.az).deg*u.deg - (m1_altaz.az).deg*u.deg
print(separation)

#3
#Generate random positions in the Galactic plane (between l=0 and l=360 and b=-1 and b=+1),
# then make a plot showing the position on the sky of these points in FK5.
rand_coords = SkyCoord([randint(1, 360), randint(1, 360), randint(1, 360)]*u.deg, [random()*randint(-1, 1), random()*randint(-1, 1), random()*randint(-1., 1)]*u.deg, frame='galactic')
rand_coords.transform_to(FK5())
rand_FK5 = rand_coords.transform_to(FK5())
first_coord = rand_FK5[0]
second_coord = rand_FK5[1]
third_coord = rand_FK5[2]
plt.plot(first_coord)
#4
#Make a plot showing the altitude above the horizon of the Sun as seen from Seattle
# over the course of today.
