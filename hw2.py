import numpy as np

# Challenge
#1. Examine the headers of the first and second HDU in the point source catalog and try adding new keywords to them or changing them
#2. Make a histogram of the G-band magnitude of the sources in the catalog - can you figure out what the upper limit for Gmag was when the table was selected?
#3. Make a plot of the position of the sources on the sky in the point source catalog, in Galactic coordinates l and b (you will need to use what we learned about SkyCoord in the previous tutorial)
#4. Try and produce a new FITS file that contains the image in the primary HDU and the above table in the second HDU


# Challenge
# 1. Make a table that contains three columns: spectral type, temperature, and radius, and incude 5 rows with fake data (or real data if you like, for example from http://www.atlasoftheuniverse.com/startype.html. Try including units on the columns that can have them.
# 2. Find the mean temperature and the maximum radius
# 3.Try and find out how to add and remove rows
# 4. Add a new column which gives the luminosity (using L=4*pi R^2 *sigma T^4)


# Challenge
# Starting from the obs table:

# 1. Make a new table that shows every other row, starting with the second row? (that is, the second, fourth, sixth, etc. rows).
# 2. Make a new table the only contains rows where name is M31


# Challenge
# Using the Mass2 table above:

# 1. Make a plot that shows j_m-h_m on the x-axis, and h_m-k_m on the y-axis

# 2. Make a new table that contains the subset of rows where the j_snr, h_snr, and k_snr columns, which give the signal-to-noise-ratio in the J, H, and K band, are greater than 10, and try and show these points in red in the plot you just made.

# 3. Make a new table (based on the full table) that contains only the RA, Dec, and the j_m, h_m and k_m columns, then try and write out this catalog into a format that you can read into another software package. For example, try and write out the catalog into CSV format, then read it into a spreadsheet software package (e.g. Excel, Google Docs, Numbers, OpenOffice). You may run into an issue at this point - if so, take a look at https://github.com/astropy/astropy/issues/7357 to see how to fix it.