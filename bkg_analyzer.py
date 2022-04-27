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
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

######################
#OPTIONS
savefiles = False
plot_profile = False
plot_spec = False
show_ima = False

#PARAMS

######################

#browse all the *.fc.fits files in a directory and its subdirectories
main_path = './Asiago_nightsky/2011/'
main_path = './'
file_ls = glob.glob(main_path+'/**/*.fc.bkg.fits', recursive= True)
names = [os.path.basename(x) for x in file_ls]

LAMBDA_bins = np.arange(3500, 8500, 100)

df = np.zeros((len(names),len(LAMBDA_bins[:-1])))


#process all the files found
file_id = 0
for name,file in zip(names,file_ls):

    #load the frame
    hdul = fits.open(file)
    hdr, data = hdul[0].header, hdul[0].data

    #take wavelength info from the hdr
    NAXIS1, NAXIS2 = hdr['NAXIS1'], hdr['NAXIS2']
    LAMBDA0, DELTA = hdr['CRVAL1'], hdr['CDELT1']
    LAMBDA_lim = hdr['UVLIM']

    #the (eventually) UV-limited wavelengths array
    LAMBDA_start = max(LAMBDA_lim, LAMBDA0)
    LAMBDA = np.arange(LAMBDA_start, LAMBDA_start+NAXIS1*DELTA, DELTA)
    if len(LAMBDA) == NAXIS1+1:
        LAMBDA = LAMBDA[:-1]

    spec = np.nanmean(data, axis=0)
    #bin_indices = np.digitize(LAMBDA, LAMBDA_bins)
    #print(bin_indices)

    hist,_ = np.histogram(LAMBDA, bins=LAMBDA_bins, weights=spec)
    df[file_id] = hist
    file_id +=1
df[29,:]=0
plt.imshow(df)
plt.xticks([0,len(LAMBDA_bins)-2],[LAMBDA_bins[0],LAMBDA_bins[-2]])
plt.yticks([0, file_id-1],[2006,2021])
plt.xlabel('binned bkg spectrum')
plt.ylabel('frame (chornological order, yr)')
plt.show()

