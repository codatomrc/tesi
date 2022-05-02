import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.signal import find_peaks, peak_widths
from datetime import datetime
import glob
import os
from wotan import flatten
from scipy.optimize import curve_fit

######################
#OPTIONS
save_plots = True
save_FITS = True
plot_profile = True
plot_spec = False
show_ima = False

#PARAMS
peak_height = 0.05 #height abouve the bkg level
data_col_frac = .75 #minimum fraction of valid pixels in a column
width_mult = 2 # interval to exclude around a source, wrt center in FWHM units
cr_width = 2.5 # cr trace spatial width
cr_prominence = 5 #treshold height wrt average column level
cr_pad = 1 # number of px to exclude around cr, in a fixed column
LAMBDA_lim = 3500 #A, limit blue wavelength
######################

#browse all the *.fc.fits files in a directory and its subdirectories
main_path = './Asiago_nightsky/2020/'
main_path = './'
file_ls = glob.glob(main_path+'/**/*.fc.fits', recursive= True)
names = [os.path.basename(x) for x in file_ls]

#initialize the .log file
f = open("bkg_extr.log", "w")
f.write('Running bkg_extr.py at '+
        datetime.now().strftime("%H:%M:%S, %Y-%m-%d")+'\n')
f.close()

#append to the log
warnings_count = 0
def log_append(message):
    f = open("bkg_extr.log", "a")
    f.write(message+'\n')
    f.close()

######################
######################

#process all the files found
for name,file in zip(names,file_ls):

    print('processing file '+name+'\n', end='\r')

    #open a FITS file
    hdul = fits.open(file)
    hdr = hdul[0].header

    #extract wavelenght information from the header
    NAXIS1, NAXIS2 = hdr['NAXIS1'], hdr['NAXIS2']
    LAMBDA0, DELTA = hdr['CRVAL1'], hdr['CDELT1']

    #generate the lambdas array
    if hdr['CTYPE1'] != 'LINEAR':
        log_append('WARNING: no linear wavelength calibration')    
    LAMBDA = np.arange(LAMBDA0, LAMBDA0+NAXIS1*DELTA, DELTA)
    if len(LAMBDA) == NAXIS1+1:
        LAMBDA = LAMBDA[:-1]

    #remove extreme blue wavelengths
    LAMBDA_start_id = len(LAMBDA)-len(LAMBDA[ LAMBDA>LAMBDA_lim])
    LAMBDA = LAMBDA[LAMBDA_start_id:]
        
    year = hdr['DATE-OBS'][:4]

    #aperture information from the hdr
    SLIT = hdr['SLIT'] #microns
    try:
        BINX, BINY = hdr['BINX'], hdr['BINY']
        TELSCALE = hdr['TELSCALE'] #arcsec/mm
        CCDSCALE = hdr['CCDSCALE'] #arcsec/px
    except KeyError:
        BINX, BINY = hdr['HBIN'], hdr['VBIN']

        log_append(' WARNING: no scale info in the hdr (using defauls)')
        
        TELSCALE = 10.70 #arcsec/mm #TO BE CHECKED!!!
        CCDSCALE = 0.60 #arcsec/px #TO BE CHECKED!!!

    SLIT_angular = SLIT/1000 * TELSCALE #slit size in arcsec
    SLIT_px = SLIT_angular / CCDSCALE / BINX #slit size in px

######################
    #bkg level estiamtion
    raw_data = hdul[0].data[:,LAMBDA_start_id:]
    raw_integr = np.sum(raw_data, axis = 1)
    x = np.arange(len(raw_integr))
        
    bkg_est = np.nanmedian(raw_integr)

######################
    #remove cosmic rays and UV noise
    data = np.copy(raw_data)

    cr_col_frac = np.zeros(len(LAMBDA)) #fraction of remaining px
    for cr_col,col in enumerate(data.T):
        col_avg = np.nanmean(data[:,cr_col])
        cr_line,_ = find_peaks(col,
                               prominence = cr_prominence*col_avg,
                               width = (0,cr_width))

        cr_widths = peak_widths(col, cr_line, rel_height=0.5)[0]

        #set left and right boundaries of the source region along the slit
        left_width = cr_line-cr_widths - cr_pad
        right_width = cr_line+cr_widths + cr_pad

        #scan each column and remove peaks
        cr_sel = np.zeros(np.shape(col), dtype=bool)
        for i in range(np.shape(col)[0]):
            for peak,width in zip(cr_line,cr_widths):
                if abs(i-peak) < width+cr_pad:
                    cr_sel[i] = True

        #counts how many pixels are left in a column
        saved_px = (NAXIS2 - np.sum(cr_sel))/NAXIS2
        cr_col_frac[cr_col] = saved_px
        if saved_px >= data_col_frac: #if enough, take the masked column
            data[cr_sel, cr_col] = np.nan
        else: #else discart the entire column
            data[:, cr_col] = 0.

######################
    #use noise/bkg info to find peaks
    integr = np.nansum(data, axis = 1)
    bkg_est = np.median(integr)
    
    #detrend: global trend (including peaks)
    _,trend_raw = flatten (x,
                            integr ,
                            method ='biweight',
                            window_length =200 ,
                            cval = 10, return_trend = True )

    #trim removing peaks, i.e. data fare above the global trend
    integr_trim = np.where(integr <= trend_raw+0.05*bkg_est,
                           integr, trend_raw)

    #detrend the trimmed data, much less sensitive to the peaks
    _,trend = flatten (x,
                       integr_trim ,
                       method ='biweight',
                       window_length =50 ,
                       cval = 10, return_trend = True )
    '''
    more plot about detrending
    plt.plot(x,integr, label='original profile')
    plt.plot(x,trend_raw, label='raw de-trend')
    plt.plot(x,integr_trim, label='trimmed profile')
    plt.plot(x,trend, label='final bkg trend estimation')
    plt.legend()
    plt.show()
    '''

    #detrend residuals: original peaks are highlighted wrt the bkg profile
    diff = integr-trend
    
    #find peaks
    peaks,properties = find_peaks(diff, height=0.05*bkg_est, width = cr_width)
    peak_FWHM = peak_widths(integr, peaks, rel_height=.5)[0]/2.

    if len(peak_FWHM)== 0:
        no_source = " WARNING: no sources were detected"
        log_append(no_source)

    #generate a boolena mask True outside the peaks
    bkg_sel = np.full(np.shape(x), True)
    for i,peak in enumerate(peaks):
        width = (int(peak_FWHM[i])+1)*width_mult
        for w in range(-width,width):
            bkg_sel[peak+w]=False
        w = width -1
        while integr[peak+w] >= trend[peak+w]:
            bkg_sel[peak+w] = False
            w += 1
        w = width
        while integr[peak-w] >=trend[peak-w]:
            bkg_sel[peak-w]=False
            w += 1
            
######################

    #plot the luminosity profile, show source and bkg regions
    if 1 == plot_profile:
        plt.title(year+'/'+name[:-8]+': wavelenght integration')
        plt.plot(raw_integr, alpha=0.2, ls='dashed', c='C1')
        plt.plot(integr, alpha=0.4) #integrated flux
        plt.scatter(x[bkg_sel],
                    integr[bkg_sel],
                    s=0.2, c='green') #select bkg

        #esimate the bkg of the filtered regions only
        bkg_est_filt = np.mean(integr[bkg_sel])

        plt.plot(x, trend_raw+0.05*bkg_est, ls='dashed', c='grey', alpha=0.5)
        ima = np.zeros(np.shape(bkg_sel))
        plt.fill_between(x, ~bkg_sel*1.1*max(integr), color='red', alpha=.1)
        
        #plt settings
        plt.ylim(min(integr)*0.9, max(integr)*1.1)
        plt.legend(['raw signal','cleaned signal',
                    'bkg signal only', 'detrended treshold',
                    f'peak regions ({width_mult}xFWHM)'])
        
        if save_plots is True:
            plt.savefig('./plots/integr/'+year+'_'+name[:-8]+'.png')
            plt.close()
        else:
            plt.show()

######################
    #integrated spectrum (along the slit)
    total = np.nanmean(data, axis = 0) #integration along the slit
    sky = np.nanmean(data[bkg_sel,:], axis = 0) #integration of bkg rows only

    total[total == 0] = np.nan

    #plot the spatially integrated spectrum of the bkg
    if 1 == plot_spec:
        plt.title(year+'/'+name[:-8]+': bkg spectrum')
        plt.plot(LAMBDA, total, label='full frame', color='gray', alpha=0.3)
        plt.plot(LAMBDA, sky, label='sky only')
        plt.legend()
        if save_plots is True:
            plt.savefig('./plots/sky_spec/'+year+'_'+name[:-8]+'.png')
            plt.close()
        else:
            plt.show()

######################
    #extract only the bkg rows
        
    ma_data = data #set masked data
    for i,row in enumerate(bkg_sel):
        #cancel data from the source rows
        if row == 0:
            ma_data[i,:] = np.nan 

    #plot as image the bkr rows only
    if (1 == show_ima) and (save_plots ==0):
        plt.title('image (sky selection only)')
        plt.imshow(ma_data, extent = [LAMBDA[0], LAMBDA[-1], NAXIS2, 0])
        plt.show()

######################
    #save masked data in a new FITS file
    if 1 == save_FITS:
        now = datetime.now()
        now_str = now.strftime("%Y-%m-%d %H:%M:%S")
        
        hdr.set('BKGEXTR', now_str, 'Time of bkg extraction')
        hdr.set('UVLIM', LAMBDA_lim, 'A')
        hdr['NAXIS1']=len(data[0])
        new_hdu = fits.PrimaryHDU(ma_data)
        new_hdul = fits.HDUList([new_hdu])
        new_hdul[0].header = hdr

        file_new = file[:-5]+'.bkg.fits'
        new_hdul.writeto(file_new, overwrite=True)

f = open("bkg_extr.log", 'r')
warnings_count = len(f.readlines())-1
if warnings_count != 0:
    print(f'WARNING: {warnings_count} warnings occurred (see the log)')
