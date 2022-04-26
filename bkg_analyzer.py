import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.coordinates as coord
from astropy import units as u
from scipy.signal import find_peaks, peak_widths, detrend
from scipy import ndimage
from datetime import datetime
import glob
import os
from wotan import flatten

######################
#OPTIONS
savefiles = False
plot_profile = False
plot_spec = False
show_ima = False

#PARAMS

######################

#browse all the *.fc.fits files in a directory and its subdirectories
main_path = './Asiago_nightsky/2006/'
#main_path = './'
file_ls = glob.glob(main_path+'/**/*.fc.bkg.fits', recursive= True)
names = [os.path.basename(x) for x in file_ls]

#process all the files found
for name,file in zip(names,file_ls):

    #load the frame
    hdul = fits.open(file)
    hdr, data = hdul[0].header, hdul[0].data

    #take wavelength info from the hdr
    NAXIS1, NAXIS2 = hdr['NAXIS1'], hdr['NAXIS2']
    LAMBDA0, DELTA = hdr['CRVAL1'], hdr['CDELT1']
    LAMBDA_lim = hdr['UVLIM']

    #the (eventually) UV-limited wavelengths array
    LAMBDA = np.arange(LAMBDA_lim, max(LAMBDA_lim, LAMBDA0)+NAXIS1*DELTA, DELTA)
    
    '''
    HERE IT GOES THE CODE TO ANALYZE THE FEATURES OF BKG SPECTRA
    '''
