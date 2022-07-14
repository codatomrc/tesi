import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, peak_widths, savgol_filter
from scipy.interpolate import UnivariateSpline, interp1d
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

plot_fits = False
save_fits = False

save_cont = False

FITS_lines = False
#PARAMS
line_res = 3 #x delta lambda, min distance to consider lines as unresolved
JD0 = 2450000
line_window = 10 #DELTA units when finding lines to be masked to estimate continuum
######################

# import lines table
lines_raw = np.genfromtxt('lines.txt', usecols=0)
ranges = np.genfromtxt('ranges.txt')

line_diff = np.diff(lines_raw)

JDs = []

#browse all the *.fc.fits files in a directory and its subdirectories
main_path = './Asiago_nightsky/2020/'
#main_path = './'
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
    close_lines = np.where(line_diff < line_res*DELTA, False, True)
    close_lines = np.insert(close_lines, 0, True)
    lines = lines_raw[close_lines]

    filename = './plots/widths/'+year+'_'+name[:-13]+'.l.txt'
    f = open(filename, 'w') if FITS_lines else 0
    f.write(f"#line\t EW") if FITS_lines else 0

    ####################

    cont_sample = np.zeros(len(LAMBDA))
    for i,Lambda in enumerate(LAMBDA):
        for Range in ranges:
            if  Lambda >= Range[0] and Lambda <= Range[1]:
                cont_sample[i] = spec[i]
    w = cont_sample != 0
    cont_sample[~w] = np.nan
    

    _,smooth_cont = flatten(LAMBDA[w], cont_sample[w], method='biweight',
                            window_length = 50, cval=1, return_trend=True)

    interp = interp1d(LAMBDA[w], smooth_cont,
                 kind = 'quadratic', fill_value="extrapolate")    
    final_cont = interp(LAMBDA)
    
    ww = np.diff(LAMBDA[w], append = 0) >= 2*DELTA   
    pippo = np.zeros(len(LAMBDA))
    pippo[~w] = np.nan
    pippo[w] = smooth_cont
    
    plt.plot(LAMBDA, spec)
    plt.plot(LAMBDA, pippo, lw=5)
    plt.plot(LAMBDA[:1500], final_cont[:1500])
    plt.show()
    

    ###################


    '''
    CONTINUUM ESTIMATION
    '''
    
    cont = np.copy(spec)
    #mask the lines from the spectrum
    widths = []
    for line in lines:
        i = int((line-LAMBDA[0])/DELTA)
        x = spec[i-line_window:i+line_window] #find peak region
        peak,prop = find_peaks(x, prominence = np.std(x)) #find lines
        width = peak_widths(x, peak, rel_height =1) #compute widths

        #mask a region of 3 x widths around the line
        try:
             peak_width = int(peak_widths(x, peak, rel_height =1)[0][0]*1.5)
        except IndexError:
            peak_width = 0 #if no peaks (i.e. faint line) mask the center only
        widths.append(peak_width)
            
        cont[i-peak_width:i+peak_width] = np.nan #perform masking

        plt.axvline(x=line, ls='--', c='k', alpha=.1)

    #flatten remaining smooth_cont (to remove noise+small lines)
    w = ~ np.isnan(cont) #take unmasked data
    _,smooth_cont = flatten(LAMBDA[w], cont[w], method='biweight',
                         window_length = 20, cval=3, return_trend=True)

    #interpolate over masked regions
    ww = ~ np.isnan(smooth_cont)
    LLAMBDA = LAMBDA[w]    
    interp = interp1d(LLAMBDA[ww], smooth_cont[ww],
                 kind = 'quadratic', fill_value="extrapolate")    
    final_cont = interp(LAMBDA)

    if save_cont is True:
        plt.plot(LAMBDA, spec, label='original')
        plt.plot(LAMBDA, cont, label='masked data')
        plt.plot(LAMBDA, final_cont, label='continuum est.')
        plt.xlabel('wavelenght [A]')
        plt.ylabel('flux [erg/cm2/s/A]')
        plt.legend()
        plt.show()
        #plt.savefig('./plots/continuum/'+year+'_'+name[:-8]+'.png', dpi=500)
        plt.close()

    '''
    LINE FIT
    '''

    #line fit and EW computation
    u_flux = u.erg / (u.cm ** 2 * u.s * u.AA) #flux units
    A = u.AA #angstrom units
    spectrum = Spectrum1D(flux=spec*u_flux, spectral_axis=LAMBDA*A)
    EWs = []

    for line,width in zip(lines,widths):
        line_init = models.Gaussian1D(amplitude=0.5*max(spec)*u_flux,
                                    mean=line*A,
                                    stddev=5.*A)
        
        line_fit = fit_lines(spectrum-final_cont, line_init,
                             window=width*DELTA*A)
        y_fit = line_fit(LAMBDA*A)

        plt.plot(LAMBDA, y_fit+final_cont*u_flux,
                 lw=0.4, ls = '-', c='C1') if plot_fits else 0

        EW = np.sum(y_fit/(final_cont*u_flux))
        EWs.append(EW)
   
    # for groups of lines the same EW is given to all the components
    EW_array = np.zeros(np.shape(lines_raw))
    EW_array[~close_lines] = np.nan #grouped lines
    EW_array[close_lines] = EWs
    for i,EW in enumerate(EW_array):
        if np.isnan(EW):
            EW_array[i] = EW_array[i-1]

        f.writelines(f"\n{lines_raw[i]}\t {EW_array[i]}") if FITS_lines else 0

    if plot_fits is True:
        plt.plot(LAMBDA, spec, lw=0.2, ls='-.')
        plt.xlabel('wavelenght [A]')
        plt.ylabel('flux [erg/cm2/s/A]')
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
            [fits.Column(name = 'line', array = lines_raw, format = 'E'),
             fits.Column(name = 'EWs', array = EW_array, format = 'E'),
             fits.Column(name = 'IsIsolated', array = close_lines, format = 'L' )])

        now = datetime.now()
        now_str = now.strftime("%Y-%m-%d %H:%M:%S")
        
        hdul.append(table_hdu)
        t_hdr = hdul[1].header
        t_hdr.set('UNITS', 'Angstrom')
        t_hdr.set('EWTIME', now_str, 'Time of EW computation')
        t_hdr.set('LINERES', now_str, 'Min dist btw lines, in DELTA units')

        #save the contimuum in a new partition too
        table_hdu = fits.BinTableHDU.from_columns(
            [fits.Column(name = 'LAMBDA', array = LAMBDA, format = 'E'),
             fits.Column(name = 'flux', array = final_cont, format = 'E')])
        hdul.append(table_hdu)

        #save the .FITS file
        file_new = file[:-12]+'.l.fits'
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

