import numpy as np
import firtez_dz as frz
import sys
from astropy.io import fits

filepath = sys.argv[1]

atm=frz.read_model(filepath)

tau_common = np.linspace(-3,1,41)
nx = atm.tem.shape[0]
ny = atm.tem.shape[1]
nz = tau_common.shape[0]
T_interp= P_interp= Bx_interp= By_interp= Bz_interp= Vx_interp= Vy_interp= Vz_interp= z_interp= tau= np.zeros((nx,ny,nz))

for i in range(nx):
    for j in range(ny):
        T_interp[i,j,:] = np.interp(tau_common, atm.tau[i,j,::-1], atm.tem[i,j,::-1])
        P_interp[i,j,:] = np.interp(tau_common, atm.tau[i,j,::-1], atm.pg[i,j,::-1])
        Bx_interp[i,j,:] = np.interp(tau_common, atm.tau[i,j,::-1], atm.bx[i,j,::-1])
        By_interp[i,j,:] = np.interp(tau_common, atm.tau[i,j,::-1], atm.by[i,j,::-1])
        Bz_interp[i,j,:] = np.interp(tau_common, atm.tau[i,j,::-1], atm.bz[i,j,::-1])
        Vx_interp[i,j,:] = np.interp(tau_common, atm.tau[i,j,::-1], atm.vx[i,j,::-1])
        Vy_interp[i,j,:] = np.interp(tau_common, atm.tau[i,j,::-1], atm.vy[i,j,::-1])
        Vz_interp[i,j,:] = np.interp(tau_common, atm.tau[i,j,::-1], atm.vz[i,j,::-1])
        z_interp[i,j,:] = np.interp(tau_common, atm.tau[i,j,::-1], atm.z[i,j,::-1])
atm = None

'''tau[:,:,:] = tau_common[None,None,:]
atm = frz.atm_model3D(nx,ny,nz)
atm.set_tem(T_interp)
atm.set_pg(P_interp)
atm.set_bx(Bx_interp)
atm.set_by(By_interp)
atm.set_bz(Bz_interp)
atm.set_vx(Vx_interp)
atm.set_vy(Vy_interp)
atm.set_vz(Vz_interp)
atm.set_tau(tau)
atm.set_z(z_interp)
atm.write_model('/home/xenosh/codes/flux_exp/firtez_test/atm_1024x1024x41_tau.bin')
'''

hdu1 = fits.PrimaryHDU(T_interp)
hdu2 = fits.ImageHDU(P_interp)
hdu3 = fits.ImageHDU(Bx_interp)
hdu4 = fits.ImageHDU(By_interp)
hdu5 = fits.ImageHDU(Bz_interp)
hdu6 = fits.ImageHDU(Vx_interp)
hdu7 = fits.ImageHDU(Vy_interp)
hdu8 = fits.ImageHDU(Vz_interp)
hdu9 = fits.ImageHDU(z_interp)
hdulist = fits.HDUList([hdu1,hdu2,hdu3,hdu4,hdu5,hdu6,hdu7,hdu8,hdu9])
hdulist.writeto(filepath[:-4]+'_remap.fits',overwrite=True)