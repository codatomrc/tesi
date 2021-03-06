import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, peak_widths
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun, get_moon, angular_separation
from astropy import units as u
from astropy.modeling import models
#from scipy.signal import find_peaks, peak_widths, detrend
from scipy import ndimage
from scipy.optimize import curve_fit
from datetime import datetime
import glob
import os
from wotan import flatten
from scipy.special import wofz
from specutils.fitting import fit_lines
from specutils.spectra import Spectrum1D
######################
#OPTIONS
plot_ranges = True
save_ranges = False

plot_fits = False
save_fits = False

#PARAMS
bin_min, bin_max = 3500, 8000
bin_size = 50
fit_window = 30 #angstrom

JD0 = 2450000
######################
Asiago = EarthLocation(lat=45.8664818*u.deg,
                       lon=11.5264273*u.deg,
                       height=1044.2*u.m)

# import lines table
lines_raw = np.genfromtxt('lines.txt', usecols=0)
line_diff = np.diff(lines_raw)
widths = np.zeros(len(lines_raw))
EWs = []
JDs = []

range_start, range_end = np.genfromtxt('ranges.txt').T
ranges = np.array([range_start,range_end]).T

#browse all the *.fc.fits files in a directory and its subdirectories
main_path = './Asiago_nightsky/2021/'
main_path = './'
file_ls = glob.glob(main_path+'/**/*.fc.bkg.fits', recursive= True)
names = [os.path.basename(x) for x in file_ls]

LAMBDA_bins = np.arange(bin_min, bin_max, bin_size)

df = np.zeros((len(names),len(LAMBDA_bins[:-1])))
airmasses, azimuths, sun_heights = [], [], []
moon_phases, moon_dist = [],[]

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
    ####

    #remove blended lines, i.e. to be considered as a single feature
    close_lines = np.where(line_diff < 5*DELTA, False, True)
    close_lines = np.insert(close_lines, 0, True)
    lines = lines_raw[close_lines]



    '''
    LINE FIT
    '''
    lvl_est = np.median(spec)
    _,trend_raw = flatten (LAMBDA,
                            spec ,
                            method ='biweight',
                            window_length =10 ,
                            cval = 2, return_trend = True )

    #trim removing peaks, i.e. data far above the global trend
    spec_trim = np.where(spec <= trend_raw, spec, trend_raw)

    #detrend the trimmed data, much less sensitive to the peaks
    _,trend = flatten (LAMBDA,
                       spec_trim ,
                       method ='biweight',
                       window_length =200 ,
                       cval = 3, return_trend = True )

    #plt.plot(LAMBDA, spec)
    #plt.plot(LAMBDA, trend)
    #plt.show()

    ######################################
    ##fino a qui sopra ok
    ######################################


    spectrum = Spectrum1D(flux=spec*u.Jy, spectral_axis=LAMBDA*u.AA)
    EWs = []
    for line in lines:
        line_init = models.Gaussian1D(amplitude=0.5*max(spec)*u.Jy,
                                    mean=line*u.AA,
                                    stddev=5.*u.AA)
        
        line_fit = fit_lines(spectrum-trend, line_init)
        y_fit = line_fit(LAMBDA*u.AA)

        EW = np.sum(y_fit/(trend*u.Jy))
        EWs.append(EW)
        
        #plt.plot(LAMBDA,y_fit+trend*u.Jy)
    plt.plot(EWs, '-o')
    

    
    '''

    #plot the regions to be fitted
    if 1==plot_ranges:
        plt.close()
        plt.plot(LAMBDA, spec, lw=0.5)
        #plt.plot(LAMBDA, fit_data*spec, lw=.5)
        for line in fit_lines:
            plt.axvline(x=line, c='C2', alpha=.2)
        plt.legend(['spectrum','fitted regions','known lines'])

        if save_ranges is True:
            plt.savefig('./plots/fit_regions/'+year+'_'+name[:-8]+'.png', dpi=500)
        else:
            plt.show()
    '''

plt.show()
#widths = widths[1:]
#cm = plt.cm.rainbow(np.linspace(0, 1, len(file_ls)))
#for i in range(len(file_ls)):
#    plt.plot(widths[i].T, color=cm[i])
plt.xlabel('line (from bluer to redder)')
plt.ylabel('EW (A)')
plt.ylim(0,+50)
sm = plt.cm.ScalarMappable(cmap='rainbow')
cbar = plt.colorbar(sm, ticks=[0,1])
cbar.set_ticklabels([2006,2021])
plt.xticks([])
plt.show()

'''

'''
#bin_indices = np.digitize(LAMBDA, LAMBDA_bins)
#print(bin_indices)

hist,_ = np.histogram(LAMBDA, bins=LAMBDA_bins, weights=spec)
df[file_id] = hist
file_id +=1

#extract other params of the frames
airmasses.append(hdr['AIRMASS'])
try:
    azimuths.append(hdr['AZIMUTH'])
except KeyError:
    azimuths.append(np.nan)

time = Time(hdr['DATE-OBS']) #utc
RA_DEC_coord = SkyCoord(hdr['RA'], hdr['DEC'], unit=(u.hourangle,u.degree))
ALT_AZ_coord = RA_DEC_coord.transform_to(AltAz(obstime=time, location=Asiago))

sun_pos = get_sun(time)
h_sun = sun_pos.transform_to(AltAz(obstime=time, location=Asiago)).alt
sun_heights.append(h_sun.value)

moon_pos = get_moon(time, location=Asiago)
phase = angular_separation(sun_pos.ra,sun_pos.dec,
                           moon_pos.ra, moon_pos.dec)/np.pi
moon_phases.append(phase.value)

moon_sep = angular_separation(RA_DEC_coord.ra,RA_DEC_coord.dec,
                           moon_pos.ra, moon_pos.dec)*180./np.pi
moon_dist.append(moon_sep.value)
'''
        

#plt.plot(LAMBDA, spec, label='spectrum')
#plt.plot(LAMBDA[, fit_data*spec, label='fitted regions')
#for line in fit_lines:
#    plt.axvline(x=line, c='C1', alpha=.2)
#plt.show()
'''
plt.plot(LAMBDA_bins[:-1], hist/max(hist))
plt.plot(LAMBDA, spec/max(spec))
plt.show()


   
#df[29,:]=0
plt.imshow(df)
plt.xticks([0,len(LAMBDA_bins)-2],[LAMBDA_bins[0],LAMBDA_bins[-2]])
plt.yticks([0, file_id-1],[2006,2021])
plt.xlabel('binned bkg spectrum')
plt.ylabel('frame (chornological order, yr)')
plt.show()


flux = np.sum(df, axis=1)
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
ax.scatter(azimuths, np.asarray(airmasses)-1., c=flux)
ax.grid(True)
ax.set_xticks([0,np.pi/2,np.pi,np.pi*3/2.], ['S', 'W', 'N', 'E'])
ax.set_xlabel('azimuth')
ax.text(0,.25,'airmass-1', rotation=0)
plt.show()

fig,ax = plt.subplots(2,2)
fig.suptitle('flux vs moon and sun positions')
ax[0,1].scatter(flux, sun_heights)
ax[0,1].set_ylabel('sun height [deg]')
ax[0,1].set_xlabel('flux [e/s/cm2]')

ax[1,0].scatter(moon_dist, flux)
ax[1,0].set_xlabel('moon sep. [deg]')
ax[1,0].set_ylabel('flux [e/s/cm2]')

ax[0,0].scatter(moon_dist, sun_heights, c=flux)
ax[0,0].set_xlabel('moon sep. [deg]')
ax[0,0].set_ylabel('sun height [deg]')

ax[1,1].scatter(moon_phases, flux)
ax[1,1].set_xlabel('moon phase')
ax[1,1].set_ylabel('flux [e/s/cm2]')

plt.tight_layout()
plt.show()
df = df/df[:,58,np.newaxis]
cov = np.cov(df.T)
plt.imshow(cov, extent = [bin_min,bin_max,bin_max,bin_min])
plt.colorbar()
plt.show()
'''
plt.close()
EWs = np.reshape(EWs, (len(names), int(len(EWs)/len(names)))).T
cm = plt.cm.Spectral(np.linspace(0, 1, len(EWs)))
for i in range(len(EWs)):
    plt.plot(np.asarray(JDs)-JD0, EWs[i].T, 'o', color=cm[i])
sm = plt.cm.ScalarMappable(cmap='turbo')
cbar = plt.colorbar(sm, ticks=[0,1])
cbar.set_ticklabels(['bluer','redder'])
plt.ylim(-100,20)
plt.xlabel(f'JD-{JD0} (d)')
plt.ylabel('EW (A)')
plt.show()
'''
