import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.coordinates as coord
from astropy import units as u
from astroquery.simbad import Simbad
from scipy.signal import find_peaks, peak_widths

hdul = fits.open('Asiago_nightsky/2006/ima_015.fc.bkg.fits')
hdr = hdul[0].header

print(hdr)


LAMBDA0 = hdr['CRVAL1']
DELTA = hdr['CDELT1']
NAXIS1 = hdr['NAXIS1']
LAMBDA = np.arange(LAMBDA0, LAMBDA0+NAXIS1*DELTA, DELTA)

data = hdul[0].data

plt.imshow(data)
plt.show()
