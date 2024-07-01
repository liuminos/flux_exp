import MilneEddington
import numpy as np
import matplotlib.pyplot as plt
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
stokes.transpose(1,2,0,3)
stok = None
print("info::read. stokes shape is: ", stokes.shape)

#normalize
#mean_continuum = np.mean(stokes[:,:,0,-10])
#stokes /= mean_continuum
#print("info::stokes cube normalized to the local continuum")

#wavelength calibration
ll1 = 6173.3354 - np.arange(51) * 12 / 1000
ll2 = 6301.5012 - np.arange(101) * 20 / 1000
ll = ll1 + ll2
#ll = 6301.0 + np.arange(201)*0.01
#noise
noise_level = 1e-3
noise = np.zeros([4,201])
noise += noise_level
#noise[0] *= 10. #noise is typicaly larger for I, because of systematics - Discuss!
noise[1] /= 5. 
noise[2] /= 5. 
noise[3] /= 5. 
print("DEBUG")
#model
regions = [[ll[:],None]]
print("DEBUG")
lines = [6173,6301,6302]
print("DEBUG")
n_threads = int(sys.argv[2]) # input number of cores as an argument
print("DEBUG")
me = MilneEddington.MilneEddington(regions, lines, nthreads=n_threads) # multiple threads
print("info::MilneEddington ready ")
#i = 0 
#j = 0 
nx = stokes.shape[0]
ny = stokes.shape[1]
model_guess = np.float64([500., 0.1, 0.1, 0.0, 0.04, 100, 0.5, 0.1, 1.0])
models_guess  = me.repeat_model(model_guess, nx, ny)
to_fit = stokes[:,:,:,:]
#me inversion
model_out, syn_out, chi2 = me.invert(models_guess, to_fit, noise, nRandom = 10, nIter=30, chi2_thres=1.0, verbose=False)
#save results
hdu1 = fits.PrimaryHDU(model_out)
#hdu2 = fits.ImageHDU(chi2)
#hdulist = fits.HDUList([hdu1,hdu2])
hdu1.writeto(filepath[:-5]+'_inverted.fits',overwrite=True)