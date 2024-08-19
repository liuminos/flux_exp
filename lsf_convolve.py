import numpy as np
from astropy.io import fits
import sys
import scipy

#load the data
#filepath = sys.argv[1]

#stokes = fits.open(filepath)[0].data
#stokes = fits.open('/dat-old/xenosh/Full_simulation/stokes_full.fits')[0].data
stokes = fits.open('/dat-old/xenosh/Full_simulation/stokes_full_converted_SOT.fits')[0].data
stokes = stokes[:,:,:,:101]
print("info::read. stokes shape is: ", stokes.shape)

# Line spread function, here use Gaussian
def Gaussian(sig,x):
    G = np.exp(-x**2/(2*sig**2))/(sig*np.sqrt(2*np.pi))
    return G
x_value = np.linspace(-20,20,41)
lsf = Gaussian(1.3,x_value)
'''
nb_filters = 6
wln_filters = np.loadtxt('wavelength.txt')
filters = np.zeros((nb_filters,len(wln_filters)))
for i in range(nb_filters):
  filters[i] = np.loadtxt('filter%i.txt' % i)
integral = np.trapz(filters[0], x=wln_filters)
filters = filters / integral
lsf = filters[3,7734:8335:10]/np.sum(filters[3,7734:8335:10])
#lsf = filters[0,7036:8637:10]/np.sum(filters[0,7036:8637:10])
'''
print("lsf prepared")

#convolve with the lsf
stokes_degraded = np.zeros(stokes.shape)
for i in range(stokes.shape[0]):
        for j in range(stokes.shape[1]):
            for k in range(stokes.shape[2]):
                stokes_degraded[i,j,k,:] = scipy.ndimage.convolve(stokes[i,j,k,:], lsf, mode='nearest')
print("convolved")
hdu = fits.PrimaryHDU(stokes_degraded)
#hdu.writeto('/dat-old/xenosh/Full_simulation/stokes_PSF_degraded_f3_shift.fits', overwrite=True)
hdu.writeto('/dat-old/xenosh/Full_simulation/stokes_SOT_degraded.fits', overwrite=True)
'''
#select the six points, here we choose 6points around the line and take average
stokesd = np.zeros((stokes_degraded.shape[0],stokes_degraded.shape[1],stokes_degraded.shape[2],6))
for i in range(stokes_degraded.shape[0]):
        for j in range(stokes_degraded.shape[1]):
            for k in range(stokes_degraded.shape[2]):
                stokesd[i,j,k,0] = np.sum(stokes_degraded[i,j,k,30:37])/7.0
                stokesd[i,j,k,1] = np.sum(stokes_degraded[i,j,k,37:44])/7.0
                stokesd[i,j,k,2] = np.sum(stokes_degraded[i,j,k,44:51])/7.0
                stokesd[i,j,k,3] = np.sum(stokes_degraded[i,j,k,51:58])/7.0
                stokesd[i,j,k,4] = np.sum(stokes_degraded[i,j,k,58:65])/7.0
                stokesd[i,j,k,5] = np.sum(stokes_degraded[i,j,k,65:72])/7.0

hdu = fits.PrimaryHDU(stokesd)
hdu.writeto('/dat-old/xenosh/Full_simulation/stokesd_PSF_degraded_f3_shift_avg.fits', overwrite=True)
'''
