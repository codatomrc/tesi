import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, peak_widths
from astropy.io import fits
from astropy import units as u
from astropy.modeling import models
from datetime import datetime
import glob
import os
from wotan import flatten
from specutils.fitting import fit_lines
from specutils.spectra import Spectrum1D
######################
#OPTIONS
plot_ranges = True
save_ranges = False

plot_fits = False
save_fits = False

#PARAMS
JD0 = 2450000
######################

# import lines table
lines_raw = np.genfromtxt('lines.txt', usecols=0)
line_diff = np.diff(lines_raw)
widths = []
JDs = []

#browse all the *.fc.fits files in a directory and its subdirectories
main_path = './Asiago_nightsky/2021/'
main_path = './'
file_ls = glob.glob(main_path+'/**/*.fc.bkg.fits', recursive= True)
names = [os.path.basename(x) for x in file_ls]

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
    year = hdr['DATE-OBS'][:4]
    JDs.append(hdr['JD'])

    #the (eventually) UV-limited wavelengths array
    LAMBDA_start = max(LAMBDA_lim, LAMBDA0)
    LAMBDA = np.arange(LAMBDA_start, LAMBDA_start+NAXIS1*DELTA, DELTA)
    if len(LAMBDA) == NAXIS1+1:
        LAMBDA = LAMBDA[:-1]
    spec = np.nanmean(data, axis=0)

    #remove blended lines, i.e. to be considered as a single feature
    close_lines = np.where(line_diff < 5*DELTA, False, True)
    close_lines = np.insert(close_lines, 0, True)
    lines = lines_raw[close_lines]

    '''
    LINE FIT
    '''
    #continuum estimation
    lvl_est = np.median(spec)
    _,trend_raw = flatten (LAMBDA,
                            spec ,
                            method ='biweight',
                            window_length =200 ,
                            cval = 10, return_trend = True )

    #trim removing peaks, i.e. data far above the global trend
    spec_trim = np.where(spec <= trend_raw+0.3*np.median(spec), spec, trend_raw)

    #detrend the trimmed data, much less sensitive to the peaks
    _,trend = flatten (LAMBDA,
                       spec_trim ,
                       method ='biweight',
                       window_length =150 ,
                       cval = 5, return_trend = True )

    #line fit and EW computation
    u_flux = u.erg / (u.cm ** 2 * u.s * u.AA) #flux units
    A = u.AA #angstrom units
    spectrum = Spectrum1D(flux=spec*u_flux, spectral_axis=LAMBDA*A)
    EWs = []
    for line in lines:
        line_init = models.Gaussian1D(amplitude=0.5*max(spec)*u_flux,
                                    mean=line*A,
                                    stddev=5.*A)
        
        line_fit = fit_lines(spectrum-trend, line_init)
        y_fit = line_fit(LAMBDA*A)

        EW = np.sum(y_fit/(trend*u_flux))
        EWs.append(EW)
        
        
    plt.plot(EWs, '-o')
    widths.append(EWs)
    

plt.show()
widths = np.asarray(widths)

#remove bad fits
neg_sel = widths <= 0.
widths[neg_sel] = np.nan


cm = plt.cm.Spectral(np.linspace(1, 0, 7500-3500))
for i in range(len(widths.T)):
    plt.plot(JDs, widths.T[i], 'o', color=cm[int(lines[i])-int(3500)])
plt.xlabel(f'Epoch [JD-{JD0}d]')
plt.ylabel('EW (A)')
#plt.ylim(0,+50)
sm = plt.cm.ScalarMappable(cmap='Spectral')
cbar = plt.colorbar(sm, ticks=[0,0.5,1])
cbar.set_ticklabels(['7500 A','5500 A','3500 A'])
plt.title('Intensities of all the detected spectral lines')
plt.xticks([])
plt.show()

