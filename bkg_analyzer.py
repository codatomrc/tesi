import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, peak_widths
from scipy.interpolate import UnivariateSpline, interp1d
from scipy.special import wofz
from scipy.optimize import curve_fit

from astropy.convolution import convolve_models
from astropy.io import fits
from astropy import units as u
from astropy.modeling import models,fitting
from astropy.modeling.functional_models import Box1D
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

show_cont = True
save_cont = True

FITS_lines = True
#PARAMS
JD0 = 2450000
line_window = 10 #DELTA units when finding lines to be masked to estimate continuum
######################

# import lines table
lines_raw = np.genfromtxt('lines.txt', usecols=0)
ranges2 = np.genfromtxt('ranges2.txt')


line_diff = np.diff(lines_raw)

JDs = []

#browse all the *.fc.fits files in a directory and its subdirectories
main_path = './Asiago_nightsky/2021/'
main_path = './Asiago_nightsky/'
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
    try:
        slit = hdr['SLIT']/1000.*hdr['TELSCALE']/hdr['CCDSCALE'] #px
    except KeyError:
        slit = hdr['SLIT']/1000.*10.70/0.63 # TO BE CHECKED!!!
    

    #the (eventually) UV-limited wavelengths array
    LAMBDA_start = max(LAMBDA_lim, LAMBDA0)
    LAMBDA = np.arange(LAMBDA_start, LAMBDA_start+NAXIS1*DELTA, DELTA)
    if len(LAMBDA) == NAXIS1+1:
        LAMBDA = LAMBDA[:-1]
    spec = np.nanmean(data, axis=0)

    #remove blended lines, i.e. to be considered as a single feature
    close_lines = np.where(line_diff < slit*DELTA, False, True)
    close_lines = np.insert(close_lines, 0, True)
    lines = lines_raw[close_lines]

    filename = './plots/widths/'+year+'_'+name[:-13]+'.l.txt'
    f = open(filename, 'w') if FITS_lines else 0
    f.write(f"#line\t EW") if FITS_lines else 0

    ####################

    '''
    CONTINUUM ESTIMATION
    '''

    cont_sample = np.zeros(len(LAMBDA)) #samples array

    x,y=[],[]
    plt.plot(LAMBDA, spec, label='original') if show_cont else 0
    #iterate over the sampling intervals
    for Range in ranges2:

        #LAMBDAs in the sampling ranges
        in_range = (LAMBDA > Range[0]) & (LAMBDA < Range[1])
        if np.sum(in_range) == 0:
            print(f'WARNING: no data in the interval from {Range[0]}A!')

        plt.plot(LAMBDA[in_range],spec[in_range], c='C2',lw=3)

        #average the signal (small interval, linear approx.)
        LAMBDA_avg = np.mean(LAMBDA[in_range])
        spec_avg = np.mean(spec[in_range])

        x.append(LAMBDA_avg)
        y.append(spec_avg)


    #interplolate the continuum from the sampling intervals
    interp = interp1d(x, y,

                 kind = 'quadratic', fill_value="extrapolate")    
    final_cont = interp(LAMBDA)
    
    #plot and save the results
    if show_cont is True:
        for line in lines:
            plt.axvline(x=line, ls='--', c='gray', lw=0.2)
        plt.plot(LAMBDA[in_range],spec[in_range], c='C2',
                 label = 'sampled regions', lw=3)
        
        plt.plot(LAMBDA, final_cont, label='continuum est.')
        plt.xlabel('wavelenght [A]')
        plt.ylabel('flux [erg/cm2/s/A]')
        plt.legend()
        if save_cont is True:
            plt.savefig('./plots/continuum/'+year+'_'+name[:-8]+'.png', dpi=500)
            plt.close()
        else:
            plt.show()
    else:
        plt.close()

    '''
    LINE FIT
    '''
    #LINE FIT
    u_flux = u.erg / (u.cm ** 2 * u.s * u.AA) #flux units
    A = u.AA #angstrom units
    spectrum = Spectrum1D(flux=spec*u_flux, spectral_axis=LAMBDA*A)
    EWs = []

    #model the line spectrum as sum of all the lines
    '''
    model = models.Gaussian1D(amplitude=0.5*max(spec)*u_flux,
                                    mean=lines[0]*A,
                                    stddev=5.*A)
    '''
    '''
    model = models.Voigt1D(x_0=lines[0], amplitude_L=0.1*max(spec),
                          fwhm_L=2*DELTA, fwhm_G=2*DELTA)
    w = Box1D(amplitude=1e15*hdr['SLIT'], width=hdr['SLIT'])
    '''
    model = models.Lorentz1D(amplitude=0.5*max(spec), x_0=lines[0], fwhm=10.)
    
    for line in lines[1:]:
        V = models.Lorentz1D(amplitude=0.5*max(spec),
                             x_0=line, fwhm=10.)
        model = model + V

        plt.axvline(x=line, lw=0.4, ls='--', c='k', alpha=0.2) if plot_fits else 0
       
    fit = fitting.LevMarLSQFitter()
    y_data = spec-final_cont
    line_fit = fit(model, LAMBDA, y_data, maxiter=1000)
    y_fit = line_fit(LAMBDA)
    '''
    #########################
    ###########################

    #line profile function
    def line_profile(x, *par):
        A,x0,alpha,gamma=list(par)
        sigma = alpha / np.sqrt(2 * np.log(2))
        
        #functional form
        V = A*np.real(wofz(((x-x0)+ 1j*gamma)/sigma/np.sqrt(2)))/sigma/np.sqrt(2*np.pi)
        w = np.full(hdr['SLIT'],1/hdr['SLIT'])
        return np.convolve(V,w, mode='same')

        
    def line_fit(x, *params):
        print('pippo')
        params=np.reshape(np.asarray(params), (int(len(params)/4),4))
        
        y = 0
        j=0
        for param in params:
            #print(j) if j == 36 else 0
            j+=1
            y += line_profile(x, *param)
        return y

    par_list = []
    for line in lines:
        par_list.append([1e-10,line,.5,.5])
    popt,pcov = curve_fit(line_profile, LAMBDA, spec-final_cont, p0=[1e-15,lines[0],10,10])
    print(popt)
    '''
    

    ###########################
    ##########################

    plt.plot(LAMBDA, y_fit+final_cont,
                lw=0.4, ls = '-', c='C1') if plot_fits else 0
    

    EW = np.sum(y_fit/(final_cont))
    EWs.append(EW)
   
    # for groups of lines the same EW is given to all the components
    EW_array = np.zeros(np.shape(lines_raw))
    EW_array[~close_lines] = np.nan #grouped lines
    EW_array[close_lines] = EWs
    for i,EW in enumerate(EW_array):
        if np.isnan(EW):
            EW_array[i] = EW_array[i-1]

        f.writelines(f"\n{lines_raw[i]}\t {EW_array[i]}") if FITS_lines else 0

    #plot the spectrum and the best fit profiles
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

    #save as a new .FITS file
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
