import numpy as np
import astropy.units as u
from astropy.units import imperial as imp
import astropy.constants as const
from astropy.time import Time
from astropy.coordinates import SkyCoord

%matplotlib inline
import matplotlib.pyplot as plt
plt.rc('image', origin='lower')
plt.rc('figure', figsize=(10, 6))
from astropy.visualization import quantity_support
quantity_support()



# Challenge: Unit conversions
time = 3 * u.Myr;
distance = np.array([1.2, 2.2, 1.7]) * u.kpc;
speed = distance / time;

# 1. Convert the speed above to miles/hour
print(speed.to(imp.mi / u.hr));

# 2. Calculate whether a pint is more than half liter.
pint = 1 * imp.pint;
half_liter = 0.5 * u.liter;
print(f"A pint is more than half a liter: {pint > half_liter}\nThis is because astropy uses the US gallon definition.");

# 3. Calculate the area of a rectangle 3 km x 5 m. Show in m^2 and in yards^2.
length = 3 * u.km;
width = 5 * u.m;
area = length * width;
print(f'Area in meters^2 is {area.si}, area in yards^2 is {area.to(imp.yd**2)}');



# Challenge: Using equivalencies
# 1. Find out more about the spectral flux equivalency and convert 12mJy to
# erg/cm^2/s/Hz and to W/m^2/Hz
flux = 12 * u.mJy;
flux1 = flux.to(u.erg / u.cm**2 / u.s / u.Hz, equivalencies=u.spectral_density(3500 * u.AA));
flux2 = flux.to(u.W / u.m**2 / u.Hz, equivalencies=u.spectral_density(3500 * u.AA));
print(f"One unit: {flux1} \nAnother unit: {flux2}");



# Challenge: Plotting the solar system potential
# 1. The gravitational potential around a point source is V=-GM/r where M is the
# mass of the point and $r$ is the radius from it. Use what we've seen above to
# make a plot of the gravitational potential (in MJ/kg) in the solar
# system between 1 and 50 AU.
r = np.linspace(1,50, 1000) * u.AU;
V = ((-const.G * const.M_sun)/r).to(u.MJ / u.kg);
plt.plot(r, V);



# Challenge
# 1. Construct a time object for the current time (note that there is a shortcut
for this)
now = Time.now();

# 2. Find a way to get an ISO 8601 string for the current time, optionally with
# 6 decimal places
now_string = Time(now, format='iso', precision=6);

# 3. Find the number of minutes that have elapsed since the start of the course
start = Time('2021-03-26 19:00:00');
difference = now - start;
print(f"It has been {difference.to(u.min):.3f} since the start of the course.");



# Challenge- finish!
# 1. Find the coordinates of the Crab Nebula (M1) in ICRS coordinates, and
# convert them to Galactic Coordinates.
m1 = SkyCoord.from_name('M1');
m1.transform_to('galactic');

# 2. Find the separation on the sky between the Crab Nebula and the star
# 'T Tauri' in degrees
t_tauri = SkyCoord.from_name('T Tauri');
(m1.separation(t_tauri)).to(u.deg);

# 3. Generate random positions in the Galactic plane (between l=0 and l=360 and
# b=-1 and b=+1), then make a plot showing the position on the sky of these
# points in FK5.


# 4. Make a plot showing the altitude above the horizon of the Sun as seen from
# Seattle over the course of today.
