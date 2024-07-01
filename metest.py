import MilneEddington
import numpy as np
from astropy.io import fits
import sys
import firtez_dz as frz

#load the data
filepath = sys.argv[1]
print("info::reading the input file...")
stok = frz.read_profile(filepath)
nw = stok.stki.shape[0]
nx = stok.stki.shape[1]
ny = stok.stki.shape[2]
stokes = np.zeros((nw,nx,ny,4))
stokes[:,:,:,0] = stok.stki[:,:,:]
stokes[:,:,:,1] = stok.stkq[:,:,:]
stokes[:,:,:,2] = stok.stku[:,:,:]
stokes[:,:,:,3] = stok.stkv[:,:,:]
stokes = stokes.transpose(1,2,3,0)
stok = None
print("info::read. stokes shape is: ", stokes.shape)

#wavelength calibration
ll1 = 6173.3354 - 0.5 + np.arange(101) * 10 / 1000
ll2 = 6301.5012 - 0.5 + np.arange(201) * 10 / 1000
ll = np.append(ll1, ll2)

#noise
noise_level = 1e-3
noise = np.zeros([4,201]) # need to change
noise += noise_level
#noise[0] *= 10. #noise is typicaly larger for I, because of systematics - Discuss!
noise[1] /= 5. 
noise[2] /= 5. 
noise[3] /= 5. 
print("DEBUG")
#model
regions = [[ll[101:],None]] # select the wavelength
print("DEBUG")
lines = [6301,6302] # lines to calculate
print("DEBUG")
n_threads = int(sys.argv[2]) # input number of cores as an argument
print("DEBUG")
me = MilneEddington.MilneEddington(regions, lines, nthreads=n_threads) # multiple threads
print("info::MilneEddington ready ")
nx = stokes.shape[0]
ny = stokes.shape[1]
model_guess = np.float64([500., 0.1, 0.1, 0.0, 0.04, 100, 0.5, 0.1, 1.0])
models_guess  = me.repeat_model(model_guess, nx, ny)
to_fit = stokes[:,:,:,101:] # !
#me inversion
model_out, syn_out, chi2 = me.invert(models_guess, to_fit, noise, nRandom = 10, nIter=30, chi2_thres=1.0, verbose=False)
#save results
hdu1 = fits.PrimaryHDU(model_out)
hdu2 = fits.ImageHDU(syn_out)
hdulist = fits.HDUList([hdu1,hdu2])
hdulist.writeto(filepath[:-4]+'_inverted_6301.fits',overwrite=True)



