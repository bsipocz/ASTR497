#Notebook 4
hdulist_psc[0].header, hdulist_psc[1].header
fits.getheader('data/gaia_lmc_psc.fits')

#Notebook 5
t1c = Table()

t1c['spectral type'] = ['O', 'B', 'A', 'F', 'G']
t1c['temperature'] = [40000, 20000, 8500, 6500, 5700]
t1c['radius'] = [10, 5, 1.7, 1.3, 1]

t1c['temperature'].unit = u.K
t1c['radius'].unit = u.solRad

t1c

np.mean(t1c['temperature'])

np.max(t1c['radius'])

t1c.add_row(['K', 4500, 0.8])

t1c

t1c.remove_row(5)

t1c

from astropy import constants as const

t1c['Luminosity'] = 4 * np.pi * (t1c['radius'] ** 2) * (t1c['temperature'] ** 4) * const.sigma_sb

t1c


obs[1:9:2]

obs.group_by('name').groups[1]