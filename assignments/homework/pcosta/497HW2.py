import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import QTable
from astropy.table import Table

#Notebook 4
#1
gaiaLMC = fits.open('../../../astropy_notebooks/data/gaia_lmc_psc.fits')
firstHDUHeader = gaiaLMC[0].header
secondHDUHeader = gaiaLMC[1].header
print(firstHDUHeader)
print(secondHDUHeader)
firstHDUHeader['BIGFISH'] = 'Always.'
secondHDUHeader['TWO'] = "No more, no less."

#2
data = gaiaLMC[1].data
gBandMag = data['phot_g_mean_mag']
num_bins = 25
n, bins, patches = plt.hist(gBandMag, num_bins, facecolor='blue')
plt.show()
maxGmag = np.amax(data['phot_g_mean_mag'])
print("Upper limit for Gmag:",maxGmag)

#3
plotCoords = SkyCoord(data['ra']*u.deg, data['dec']*u.deg, frame='fk5')
plotGalactic = plotCoords.transform_to('galactic')
lCoord = plotGalactic.l
bCoord = plotGalactic.b
plt.plot(lCoord, bCoord, 'o', markersize=2)
plt.xlabel("l Coordinate (Degrees)")
plt.ylabel("b coordinate (Degrees)")
plt.show()

#4
hduImage = fits.ImageHDU(fits.open('../../../astropy_notebooks/data/LMCDensFits1k.fits'))
hduTable = fits.BinTableHDU(fits.open('../../../astropy_notebooks/data/gaia_lmc_psc.fits'))
hduImage.writeto('myFits.fits')
hduTable.writeto('myFits.fits')


#Notebook 5
#1.1
table1 = QTable()
table1['Spectral Type'] = ['O', 'B', 'A', 'F', 'G']
table1['Temperature'] = [40000, 20000,8500, 6500, 5700]
table1['Temperature'].unit = u.K
table1['Radius'] = [10, 5, 1.7, 1.3, 1.0]
table1['Radius'].unit = u.Msun
print(table1)

#1.2
meanTemp = np.mean((np.array(table1['Temperature'])))*u.K
meanRadius = np.mean(np.array(table1['Radius']))*u.Msun
print("Mean Temp:", meanTemp, "Mean Radius:", meanRadius)

#1.3
#delete a row
table1.remove_row(2)
#add a row
table1.add_row(['K', 4500*u.K, 0.8*u.Msun])
print(table1)

#1.4
stefBoltz = 5.67*10**(-8)*u.W/((u.m**2)*u.K**4)
LumData = np.array(4*np.pi*(table1['Radius']**2)*stefBoltz*(table1['Temperature']**4))
table1['Luminosity'] = LumData
print(table1)

#2.1
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
otherObs = obs[1::2]
print(otherObs)

#2.2
onlyM31 = obs[obs['name'] == 'M31']
print(onlyM31)

#3.1
t6 = Table.read('../../../astropy_notebooks/data/2mass.tbl', format='ascii.ipac')
xCoords = t6['j_m'] - t6['h_m']
yCoords = t6['h_m'] - t6['k_m']
plt.xlabel('j_m - h_m')
plt.ylabel('h_m - k_m')
plt.plot(xCoords, yCoords, 'o', markersize=2)

#3.2
#get the rows where j_snr > 10
jRows = t6[t6['j_snr'] > 10]
#get the rows where j_snr > 10 and h_snr > 10
jhRows = jRows[jRows['h_snr'] > 10]
#get the rows where j_,h_,k_snr > 10 (this can likely be consolidated into one line)
jhkRows = jhRows[jhRows['k_snr'] > 10]
jhkXCoords = jhkRows['j_m'] - jhkRows['h_m']
jhkYCoords = jhkRows['h_m'] - jhkRows['k_m']
plt.plot(jhkXCoords, jhkYCoords, 'ro', markersize=2)
plt.show()

#3.3
aNewTable = Table()
#the following 5 lines can likely be consolidated to 1, but currently unsure how
aNewTable['ra'] = t6['ra']
aNewTable['dec'] = t6['dec']
aNewTable['j_m'] = t6['j_m']
aNewTable['h_m'] = t6['h_m']
aNewTable['k_m'] = t6['k_m']
aNewTable.write('497NewTable.csv')




