import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.coordinates as coord
from astropy import units as u
from scipy.signal import find_peaks, peak_widths
from datetime import datetime
import glob

######################

#browse all the *.fc.fits files in a directory and its subdirectories
main_path = '.'
file_ls = glob.glob(main_path+'/**/*.fc.fits', recursive= True)

######################
######################

#process all the files found
for file in file_ls:

    #open a FITS file
    hdul = fits.open(file)
    hdr = hdul[0].header

    #extract wavelenght information from the header
    LAMBDA0 = hdr['CRVAL1']
    DELTA = hdr['CDELT1']
    NAXIS1 = hdr['NAXIS1']
    LAMBDA = np.arange(LAMBDA0, LAMBDA0+NAXIS1*DELTA, DELTA)
    NAXIS2 = hdr['NAXIS2']

######################
    #integrate along the wavelenghts (dispersion direction)
    integr = 0
    for i in range(2047):
        integr = integr + hdul[0].data[:,i]

    #find peaks and measure their widths
    peaks,properties = find_peaks(integr, prominence=max(integr)/200.)
    peak_FWHM = peak_widths(integr, peaks, rel_height=0.5)[0]

    #set left and right boundaries of the source region along the slit
    width_mult = 5
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
    if left_width[0] < 0:
        left_width[0] = 0
    if right_width[-1] > NAXIS2:
        right_width[-1] = NAXIS2

######################
    #mask the source/background regions
    sign_sel = np.zeros(np.shape(integr), dtype=bool)
    for i in range(len(integr)):
        for peak,width in zip(peaks,peak_FWHM):
            if abs(i-peak) < width*width_mult:
                sign_sel[i] = True
    bkg_sel = ~sign_sel

    #plot the integrated flux, show source and bkg regions
    if 1 == False:
        plt.title('wavelenght integration')
        plt.plot(integr, alpha=0.4) #integrated flux
        plt.scatter(np.arange(len(integr))[bkg_sel],
                    integr[bkg_sel],
                    s=0.2, c='green') #select bkg

        for i in range(len(left_width)): #show all the source regions
            plt.axvspan(left_width[i],
                        right_width[i],
                        alpha=0.1, color='red')

        plt.legend(['integrated signal','bkg signal only', f'peak regions ({width_mult}xFWHM)'])
        plt.show()

######################
    #integrated spectrum (along the slit)
    spectrum = hdul[0].data #original data
    total = np.sum(spectrum, axis = 0) #integration along the slit
    sky = np.sum(spectrum[bkg_sel,:], axis = 0) #integration of bkg rows only

    #plot the spatially integrated spectrum of the bkg
    if 1 == False:
        plt.title('spatially integrated background spectrum')
        plt.plot(LAMBDA, total, label='full frame', color='gray', alpha=0.3)
        plt.plot(LAMBDA, sky, label='sky only')
        plt.ylim(min(sky),1.2*max(sky))
        plt.legend()
        plt.show()

######################
    #extract only the bkg rows
        
    ma_spectrum = spectrum #set masked data
    for i,row in enumerate(bkg_sel):
        #cancel data from the source rows
        if row == 0:
            ma_spectrum[i,:] = np.nan 

    #plot as image the bkr rows only
    if 1 == False:
        plt.title('image (sky selection only)')
        plt.imshow(ma_spectrum, extent = [LAMBDA[0], LAMBDA[-1], NAXIS2, 0])
        plt.show()

######################
    #save masked data in a new FITS file
    if 1 == False:
        now = datetime.now()
        now_str = now.strftime("%Y-%m-%d %H:%M:%S")
        
        hdr.set('BKGEXTR', now_str, 'Time of bkg extraction')
        new_hdu = fits.PrimaryHDU(ma_spectrum)
        new_hdul = fits.HDUList([new_hdu])
        new_hdul[0].header = hdr

        file_new = file[:-5]+'.bkg.fits'
        new_hdul.writeto(file_new, overwrite=True)
        print(file_new,' saved')
