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
savefiles = True

######################

#wavelenght integration function
def integrate(image, LAMBDA):
    integr = 0
    if len(LAMBDA) == NAXIS1+1:
        LAMBDA = LAMBDA[:-1]
    for i in range(len(LAMBDA)):
        integr = integr + image[:,i]
    return integr

#browse all the *.fc.fits files in a directory and its subdirectories
#main_path = './Asiago_nightsky/2006/'
main_path = './'
file_ls = glob.glob(main_path+'/**/*.fc.fits', recursive= True)
names = [os.path.basename(x) for x in file_ls]

######################
######################

#process all the files found
for name,file in zip(names,file_ls):

    print('processing file '+name+'\n', end="\r")

    #open a FITS file
    hdul = fits.open(file)
    hdr = hdul[0].header

    #extract wavelenght information from the header
    LAMBDA0 = hdr['CRVAL1']
    DELTA = hdr['CDELT1']
    NAXIS1 = hdr['NAXIS1']
    LAMBDA = np.arange(LAMBDA0, LAMBDA0+NAXIS1*DELTA, DELTA)
    NAXIS2 = hdr['NAXIS2']
    year = hdr['DATE-OBS'][:4]

    #aperture information from the hdr
    SLIT = hdr['SLIT'] #microns
    try:
        BINX = hdr['BINX']
        BINY = hdr['BINY']
        TELSCALE = hdr['TELSCALE'] #arcsec/mm
        CCDSCALE = hdr['CCDSCALE'] #arcsec/px
    except KeyError:
        BINX = hdr['HBIN']
        BINY = hdr['VBIN']
        print('Using default CCD and focal scales')
        
        TELSCALE = 10.70 #arcsec/mm #TO BE CHECKED!!!
        CCDSCALE = 0.60 #arcsec/px #TO BE CHECKED!!!

    SLIT_angular = SLIT/1000 * TELSCALE #slit size in arcsec
    SLIT_px = SLIT_angular / CCDSCALE / BINX #slit size in px

    #PSF = max(LAMBDA)/(1.22*1e10) *206265 # PSF size in arcsec
    #PSF_px = PSF / CCDSCALE / BINY #PSF size in px
    #print(PSF_px)

    #gain and ron info from the hdr
    gain = hdr['GAIN']
    ron = hdr['RDNOISE']

######################
    #bkg level estiamtion
    raw_data = hdul[0].data
    
    if len(LAMBDA) == NAXIS1+1:
        LAMBDA = LAMBDA[:-1]
    raw_integr = np.sum(raw_data, axis = 1)       
    x = np.arange(len(raw_integr))
    
    #signal detrend
    flat,trend = flatten(
        x,
        raw_integr,
        method='biweight',
        window_length=10*SLIT_px,
        cval = 1, return_trend = True)

    #esimate the bkg level and noise amplitude
    bkg_est = np.nanmean(trend)
    noise = np.nanstd(raw_integr-trend)

######################
    #remove cosmic rays and UV noise
    data = np.copy(raw_data)

    pad = 1 # number of px outside cosmic rays
    for cr_col,col in enumerate(data.T):
        cr_line,_ = find_peaks(col,
                               prominence = bkg_est/len(raw_data[1])*5,
                               width = (0,3))

        cr_widths = peak_widths(col, cr_line, rel_height=0.5)[0]

        #set left and right boundaries of the source region along the slit
        width_mult = 3
        left_width = cr_line-cr_widths - pad
        right_width = cr_line+cr_widths + pad

        cr_sel = np.zeros(np.shape(col), dtype=bool)
        for i in range(np.shape(col)[0]):
            for peak,width in zip(cr_line,cr_widths):
                if abs(i-peak) < width*width_mult:
                    cr_sel[i] = True

        data[cr_sel, cr_col] = np.nan

######################
    #use noise/bkg info to find peaks

    integr = np.nansum(data, axis = 1)
    
    peaks,properties = find_peaks(integr, prominence=noise, width = 3)
    peak_FWHM = peak_widths(integr, peaks, rel_height=0.5)[0]

    #set left and right boundaries of the source region along the slit
    width_mult = 3
    left_width = peaks-peak_FWHM*width_mult
    right_width = peaks+peak_FWHM*width_mult



######################
    #remove overlapping ranges

    #compress left and right boundaries into a single array
    widths = np.array([left_width.T, right_width.T]).T
    widths = np.reshape(widths, 2*len(left_width))

    #find and mark overlaps
    for i in range(len(widths)-1):
        if widths[i] > widths[i+1]:
            widths[i] = np.nan
            widths[i+1] = np.nan
    #shrink the array to unmasked values
    widths=widths[~np.isnan(widths)]

    #reshape the array into the original two
    left_width, right_width = np.reshape(widths, (int(len(widths)/2),2)).T

    #remove values beyond the CCD size limits
    try:
        if left_width[0] < 0:
            left_width[0] = 0
        if right_width[-1] > NAXIS2:
            right_width[-1] = NAXIS2
    except IndexError:
        print('no sources were detected!')

######################
    #mask the source/background regions
    sign_sel = np.zeros(np.shape(integr), dtype=bool)
    for i in range(len(integr)):
        for peak,width in zip(peaks,peak_FWHM):
            if abs(i-peak) < width*width_mult:
                sign_sel[i] = True
    bkg_sel = ~sign_sel

    #plot the integrated flux, show source and bkg regions
    if 1 == True:
        plt.title(year+'/'+name[:-8]+': wavelenght integration')
        plt.plot(raw_integr, alpha=0.2, ls='dashed', c='C1')
        plt.plot(integr, alpha=0.4) #integrated flux
        plt.scatter(x[bkg_sel],
                    integr[bkg_sel],
                    s=0.2, c='green') #select bkg

        #esimate the bkg of the filtered regions only
        bkg_est_filt = np.mean(integr[bkg_sel])

        plt.axhline(y=bkg_est_filt, ls='dashed', c='grey', alpha=0.5)
        for i in range(len(left_width)): #show all the source regions
            plt.axvspan(left_width[i],
                        right_width[i],
                        alpha=0.1, color='red')
        #plt settings
        plt.ylim(min(integr)*0.9, max(integr)*1.1)
        plt.legend(['raw signal','cleaned signal',
                    'bkg signal only', 'bkg level',
                    f'peak regions ({width_mult}xFWHM)'])
        
        if savefiles is True:
            plt.savefig('./plots/integr/'+year+'_'+name[:-8]+'.png')
            plt.close()
        else:
            plt.show()

######################
    #integrated spectrum (along the slit)
    total = np.nansum(data, axis = 0) #integration along the slit
    sky = np.nansum(data[bkg_sel,:], axis = 0) #integration of bkg rows only

    #plot the spatially integrated spectrum of the bkg
    if 1 == True:
        plt.title(year+'/'+name[:-8]+': bkg spectrum')
        plt.plot(LAMBDA, total, label='full frame', color='gray', alpha=0.3)
        plt.plot(LAMBDA, sky, label='sky only')
        plt.ylim(min(sky),1.2*max(sky))
        plt.legend()
        if savefiles is True:
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
    if 1 == False:
        plt.title('image (sky selection only)')
        plt.imshow(ma_data, extent = [LAMBDA[0], LAMBDA[-1], NAXIS2, 0])
        plt.show()

######################
    #save masked data in a new FITS file
    if 1 == True:
        now = datetime.now()
        now_str = now.strftime("%Y-%m-%d %H:%M:%S")
        
        hdr.set('BKGEXTR', now_str, 'Time of bkg extraction')
        new_hdu = fits.PrimaryHDU(ma_data)
        new_hdul = fits.HDUList([new_hdu])
        new_hdul[0].header = hdr

        file_new = file[:-5]+'.bkg.fits'
        new_hdul.writeto(file_new, overwrite=True)
        print(file_new,' saved')
