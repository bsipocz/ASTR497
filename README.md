#### ASTR 497  - 2021 Spring

The course focuses on tools in the Astropy ecosystem, including the core library and its
affiliated packages. It provides an introduction to the variety of functionality which is
available inside Astropy with hands-on Jupyter notebooks.

Topics include core functionality of units, quantities & constants; coordinates; FITS, ASCII,
and Astropy tables; modeling; WCS and images; timeseries. From the wider ecosystem, we will also
cover topics to access data with astroquery and to perform point-source photometry with photutils.

We also cover topics related to open source contribution workflows, setting up
repositories with testing and documentation, including hands-on exercises of opening pull
requests and performing peer code review.



#### Location
 * When: Fridays 12:00 -- 14:40 Spring quarter 2021
 * Where: Zoom




#### Requirements

To run all the notebooks on your laptop, you will need the following software
installed:

* Python 3.6 or later
* Jupyter Notebook
* Numpy
* Matplotlib
* SciPy
* Astropy 4.0 or later
* reproject
* photutils
* regions
* astropy-healpix
* astroquery
* APLpy
* PyVO


#### Assingments

Assingments are due 2 weeks after the relevant notebooks are covered in class. The exception is the first
assigment; it is due 2 weeks after the PR workflow (that is expected to be used for submission) is covered.
Early submissions are encouraged, however deadlines follow open source deadline tradition and are in AoE timezone.

* Units - Time - Coordinates: 2021-04-30
* Fits - Table - Unified I/O: 2021-05-14
* WCS - WCSAxes: 2021-05-21
* NDData - Convolution - Reproject - Statistic: 2021-05-28
* Modeling - Photometry: 2021-06-04


Form of submission:

* For each homework assigment open a pull request to this repository, with code placed in the
``assignments/homework/<yourname>`` directory.
Have the exercises submitted as ``.py`` files rather than notebooks, but further structuring is your
choice (one file per topic, or one file per assigment, etc.). Performing code review on others' submission
may be part of your next assigment.

* Please do not include any changes to the class notebooks in your submission, neither submit temporary files,
draft notebooks, etc.
