import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun, get_moon, angular_separation
from astropy import units as u
#from scipy.signal import find_peaks, peak_widths, detrend
from scipy import ndimage
from scipy.optimize import curve_fit
from datetime import datetime
import glob
import os
######################
#OPTIONS
savefiles = False
plot_profile = False
plot_spec = False
show_ima = False

#PARAMS
bin_min, bin_max = 3500, 8000
bin_size = 50
fit_window = 10 #anggrom
######################
Asiago = EarthLocation(lat=45.8664818*u.deg,
                       lon=11.5264273*u.deg,
                       height=1044.2*u.m)

# import lines table
lines = np.genfromtxt('lines.txt', usecols=0)
widths = np.zeros(len(lines))

#define gaussian function for line fit
def gauss(x, H, A, x0, sigma):
    return H + A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))

#browse all the *.fc.fits files in a directory and its subdirectories
main_path = './Asiago_nightsky/2011/'
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

    #the (eventually) UV-limited wavelengths array
    LAMBDA_start = max(LAMBDA_lim, LAMBDA0)
    LAMBDA = np.arange(LAMBDA_start, LAMBDA_start+NAXIS1*DELTA, DELTA)
    if len(LAMBDA) == NAXIS1+1:
        LAMBDA = LAMBDA[:-1]

    spec = np.nanmean(data, axis=0)

    EWs = []
    for line in lines:
        line_id = np.argmin(abs(LAMBDA-line))
        x = LAMBDA[line_id-fit_window: line_id+fit_window]
        y = spec[line_id-fit_window: line_id+fit_window]
        params,_ = curve_fit(gauss,x,y, p0=[1e-16, 1e-15, line, 5])
        y_fit =gauss(x, *params)
        #plt.plot(x,y_fit)
        #plt.axvline(x=line, c='C1', alpha=.2)

        #equivalent width, area from the analitycal gaussian integral
        EW = 1/params[0]*params[1]*np.sqrt(2*np.pi)*params[-1]
        EWs.append(EW)
    widths = np.vstack([widths, EWs])
    #plt.plot(LAMBDA, spec)
widths = widths[1:]
cm = plt.cm.rainbow(np.linspace(0, 1, len(file_ls)))
for i in range(len(file_ls)):
    plt.plot(widths[i].T, color=cm[i])
plt.xlabel('line (from bluer to redder)')
plt.ylabel('EW (A)')
plt.ylim(0,+50)
sm = plt.cm.ScalarMappable(cmap='rainbow')
cbar = plt.colorbar(sm, ticks=[0,1])
cbar.set_ticklabels([2006,2021])
plt.xticks([])
plt.show()



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
        

plt.plot(LAMBDA, spec/max(spec))
for line in lines:
    plt.axvline(x=line, c='C1', alpha=.2)
plt.show()
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
