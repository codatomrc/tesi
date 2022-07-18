import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import units as u
from astropy.table import Table
from datetime import datetime
import glob
import os
######################
#OPTIONS

plot_fits = False
save_fits = False

FITS_lines = False
#PARAMS
line_res = 2 #x delta lambda, 
JD0 = 2450000
######################

blue_cont = []
HPS_cont = []
ctr_cont = []

#browse all the *.fc.fits files in a directory and its subdirectories
main_path = './Asiago_nightsky/2006/'
main_path = './'
file_ls = glob.glob(main_path+'/**/*.l.fits', recursive= True)
names = [os.path.basename(x) for x in file_ls]

#initialize array that contain line information
data_array = np.zeros((len(file_ls)+1,
                       len(fits.open(file_ls[0])[1].data)+1))

#process all the files found
for i,file in enumerate(file_ls):

    #load the frame
    hdul = fits.open(file)
    hdr = hdul[0].header

    #take wavelength info from the hdr
    year = hdr['DATE-OBS'][:4]

    data_array[i+1,0] = hdr['JD']

    line_widths = hdul[1].data

    for j,line_tuple in enumerate(line_widths):
        data_array[i+1,j+1] = line_tuple[1]
        data_array[0,j+1] = line_tuple[0]

    continuum_data = hdul[2].data

    #HPS lamp continuum
    cont_val = 0
    for datum in continuum_data:
        if datum[0] >= 5600 and datum[0] <= 6500:
            cont_val += datum[1]
    HPS_cont.append(cont_val)

    #blue continuum
    cont_val = 0
    for datum in continuum_data:
        if datum[0] >= 4200 and datum[0] <= 4800:
            cont_val += datum[1]
    blue_cont.append(cont_val)

    #control band continuum
    cont_val = 0
    for datum in continuum_data:
        if datum[0] >= 7000 and datum[0] <= 8000:
            cont_val += datum[1]
    ctr_cont.append(cont_val)    
            

'''
PLOT INTERESTING LINES
'''
def plot_lines(a, line_list):
    x = np.zeros(len(a)-1)
    y = np.zeros((len(line_list), len(a)-1))
    for j,line in enumerate(line_list): #iterate over the lines 
        line_sel = np.round(a[0], decimals = 2) == line
        for i in range(1,len(a)): #iterate over the frames
            x[i-1] = a[i,0]
            y[j,i-1] = a[i, line_sel]
    return x,y

Hg_lines = [3650.15, 4046.56, 4077.84, 4358.34, 4670.83, 5354.03, 5460.75, 5675.81, 5769.61, 5790.67]
Na_lines = [4978.54, 5148.84, 5682.63, 6154.23] #skipped NaD # 5892
OI_lines = [5577.34, 6300.30]

plt.title('Hg lines')
x,y = plot_lines(data_array, Hg_lines)
for i,yy in enumerate(y[1:]):
    order = x.argsort()
    plt.plot(x[order]-JD0,yy[order], label=str(Hg_lines[i]))
plt.xlabel(f'epoch-{JD0} [JD]')
plt.ylabel('EW [A]')
plt.legend()
plt.show()


plt.title('Na lines')
x,y = plot_lines(data_array, Na_lines)
for i,yy in enumerate(y[1:]):
    order = x.argsort()
    plt.plot(x[order]-JD0,yy[order], label=str(Na_lines[i]))
plt.xlabel(f'epoch-{JD0} [JD]')
plt.ylabel('EW [A]')
plt.legend()
plt.show()


plt.title('Telluric lines ([0I] transitions)')
x,y = plot_lines(data_array, OI_lines)
for i,yy in enumerate(y):
    order = x.argsort()
    plt.plot(x[order]-JD0,yy[order], label=str(OI_lines[i]))
plt.xlabel(f'epoch-{JD0} [JD]')
plt.ylabel('EW [A]')
plt.legend()
plt.show()

#continuums
norm_blue = np.asarray(blue_cont)/np.asarray(ctr_cont)
norm_HPS = np.asarray(HPS_cont)/np.asarray(ctr_cont)

plt.title('LED blue continuum')
plt.plot(x[order]-JD0, norm_blue[order])
plt.xlabel(f'epoch-{JD0} [JD]')
plt.ylabel('normalized flux')
plt.show()

plt.title('HPS lamps continuum')
plt.plot(x[order]-JD0, norm_HPS[order], c='C1')
plt.xlabel(f'epoch-{JD0} [JD]')
plt.ylabel('normalized flux')
plt.show()

#ratio
plt.title('HPS/blue continuums (the higher the redder)')
plt.plot(x[order]-JD0, norm_HPS[order]/norm_blue[order])
plt.axhline(y=1, ls='--', c='gray')
plt.xlabel(f'epoch-{JD0} [JD]')
plt.ylabel('normalized flux')
plt.show()

