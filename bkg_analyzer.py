import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, peak_widths
from astropy.io import fits
from astropy import units as u
from astropy.modeling import models
from astropy.table import Table
from datetime import datetime
import glob
import os
from wotan import flatten
from specutils.fitting import fit_lines
from specutils.spectra import Spectrum1D
######################
#OPTIONS

plot_fits = True
save_fits = True

FITS_lines = False
#PARAMS
JD0 = 2450000
######################

# import lines table
lines_raw = np.genfromtxt('lines.txt', usecols=0)
line_diff = np.diff(lines_raw)
widths = []
JDs = []

#browse all the *.fc.fits files in a directory and its subdirectories
main_path = './Asiago_nightsky/2006/'
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

    filename = './plots/widths/'+year+'_'+name[:-13]+'.l.txt'
    f = open(filename, 'w') if FITS_lines else 0
    f.write(f"#line\t EW") if FITS_lines else 0

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

        plt.plot(LAMBDA, y_fit+trend*u_flux,
                 lw=2, ls = '-', c='C1') if plot_fits else 0

        EW = np.sum(y_fit/(trend*u_flux))
        EWs.append(EW)

        f.writelines(f"\n{line}\t {EW}") if FITS_lines else 0

    plt.plot(LAMBDA, spec, lw=1, ls='-.') if plot_fits else 0
    if (save_fits is True) and (plot_fits is True):
        plt.savefig('./plots/line_fit/'+year+'_'+name[:-8]+'.png', dpi=500)
    elif plot_fits is True:
        plt.show()
    plt.close()   
        
    #plt.plot(EWs, '-o')
    f.close() if FITS_lines else 0

    if FITS_lines is True:
        #save new FIT file with with EW in a partition
        table_hdu = fits.BinTableHDU.from_columns(
            [fits.Column(name = 'line', array = lines, format = 'E'),
             fits.Column(name = 'EWs', array = EWs, format = 'E')])

        now = datetime.now()
        now_str = now.strftime("%Y-%m-%d %H:%M:%S")
        
        hdul.append(table_hdu)
        t_hdr = hdul[1].header
        t_hdr.set('UNITS', 'Angstrom')
        t_hdr.set('EWTIME', now_str, 'Time of EW computation')
        
        
        x = file[:-12]+'.l.bkg.fits'
        hdul.writeto(file_new, overwrite=True)
    
    

#plt.show()
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

