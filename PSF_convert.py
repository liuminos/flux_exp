import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.convolution as apconv
import sys

#load data
filepath = sys.argv[1]

stokes = fits.open(filepath)[0].data
print("info::read. stokescube shape is: ", stokes.shape)

#nomarlize
mean_continuum = np.mean(stokes[:,:,0,-10])
stokes /= mean_continuum
print("info::stokes cube normalized to the local continuum")

#wavelength calibration
ll = np.arange(201)*0.01+6301.0

#prepare the PSF
#D = 0.20 #m
D = int(sys.argv[2])
llambda = 500e-9
pix_size = 16e3 #m
pixel_scale_arcsec = 16e3 / 725e3
diff_limit_arcsec = 1.22 * llambda / D * 206265
diff_limit_pixels = diff_limit_arcsec / pixel_scale_arcsec
PSF = apconv.AiryDisk2DKernel(diff_limit_pixels, mode='oversample')
print("info::PSF prepared")

#convolve PSF with data and bin
convolve = np.zeros([1536,1536,4,201])
cb = np.zeros([768,768,4,201])

for i in range(4):
    for j in range(201):
        convolve[:,:,i,j] = apconv.convolve_fft(stokes[:,:,i,j],PSF,boundary='wrap',normalize_kernel=True)
        cb[:,:,i,j] = np.sum(convolve[:,:,i,j].reshape(768,2,768,2),axis=(1,3))/4.0

#save
hdu1 = fits.PrimaryHDU(cb)
#hdu2 = fits.ImageHDU(convolve)
#hdulist = fits.HDUList([hdu1,hdu2])
#hdulist.writeto(filepath[:-5]+'_converted_0.1.fits', overwrite=True)
hdu1.writeto(filepath[:-5]+'_converted_'+str(D)+'.fits', overwrite=True)