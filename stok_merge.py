import numpy as np
import firtez_dz as frz
from astropy.io import fits


stok = frz.read_profile('/dat/xenosh/Full_simulation/stokes_00.bin')
nw = stok.stki.shape[0]
nx = stok.stki.shape[1]
ny = stok.stki.shape[2]
stokes1 = np.zeros((nw,nx,ny,4))
stokes1[:,:,:,0] = stok.stki[:,:,:]
stokes1[:,:,:,1] = stok.stkq[:,:,:]
stokes1[:,:,:,2] = stok.stku[:,:,:]
stokes1[:,:,:,3] = stok.stkv[:,:,:]
stokes1 = stokes1.transpose(1,2,3,0)
stok = None 
print('info::first file done',stokes1.shape)

stok = frz.read_profile('/dat/xenosh/Full_simulation/stokes_01.bin')
stokes2 = np.zeros((nw,nx,ny,4))
stokes2[:,:,:,0] = stok.stki[:,:,:]
stokes2[:,:,:,1] = stok.stkq[:,:,:]
stokes2[:,:,:,2] = stok.stku[:,:,:]
stokes2[:,:,:,3] = stok.stkv[:,:,:]
stokes2 = stokes2.transpose(1,2,3,0)
stok = None
print('info::2nd file done',stokes2.shape)
stok = frz.read_profile('/dat/xenosh/Full_simulation/stokes_10.bin')
stokes3 = np.zeros((nw,nx,ny,4))
stokes3[:,:,:,0] = stok.stki[:,:,:]
stokes3[:,:,:,1] = stok.stkq[:,:,:]
stokes3[:,:,:,2] = stok.stku[:,:,:]
stokes3[:,:,:,3] = stok.stkv[:,:,:]
stokes3 = stokes3.transpose(1,2,3,0)
stok = None
print('info::3rd file done',stokes3.shape)
stok = frz.read_profile('/dat/xenosh/Full_simulation/stokes_11.bin')
stokes4 = np.zeros((nw,nx,ny,4))
stokes4[:,:,:,0] = stok.stki[:,:,:]
stokes4[:,:,:,1] = stok.stkq[:,:,:]
stokes4[:,:,:,2] = stok.stku[:,:,:]
stokes4[:,:,:,3] = stok.stkv[:,:,:]
stokes4 = stokes4.transpose(1,2,3,0)
stok = None
print('info::4th file done',stokes4.shape)
stokes = np.append(np.append(stokes1,stokes3,axis=0),np.append(stokes2,stokes4,axis=0),axis=1)
print('info::merge file done',stokes.shape)


#save 
hdu = fits.PrimaryHDU(stokes[:,:,:,:])
hdu.writeto('/dat/xenosh/Full_simulation/stokes_full.fits',overwrite=True)