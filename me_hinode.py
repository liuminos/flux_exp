import MilneEddington
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys

#load the data
filepath = sys.argv[1]

print("info::reading the input file...")
stokes = fits.open(filepath)[0].data
print("info::read. stokescube shape is: ", stokes.shape)
#doubleGrid
nx, ny, ns, nw = stokes.shape
obs = np.zeros((nx,ny,ns,nw*2))
obs[:,:,:,0::2] = stokes
print("info::doubleGrid done.")
#normalize
#mean_continuum = np.mean(stokes[:,:,0,-10])
#stokes /= mean_continuum
#print("info::stokes cube normalized to the local continuum")

#wavelength calibration
ll = 6302.08 +np.arange(112*2)*0.0215/2-56.5*0.0215
tr = np.float64([0.00240208, 0.00390950, 0.0230995, 0.123889, 0.198799,0.116474,0.0201897,0.00704875,0.00277027]) # source A. Asensio

#noise
noise = np.zeros([4,112*2]) +1.e32
noise[:,0::2] = 1.e-3
noise[1:3] /=4.0
noise[3] /= 3.0
 
print("DEBUG")
#model
regions = [[ll[:],tr/tr.sum()]]
print("DEBUG")
lines = [6301,6302]
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
to_fit = obs[:,:,:,:]
#me inversion
model_out, syn_out, chi2 = me.invert(models_guess, to_fit, noise, nRandom = 10, nIter=30, chi2_thres=1.0, verbose=False)
#save results
hdu1 = fits.PrimaryHDU(model_out)
#hdu2 = fits.ImageHDU(chi2)
#hdulist = fits.HDUList([hdu1,hdu2])
hdu1.writeto(filepath[:-5]+'_inverted.fits',overwrite=True)