import numpy as np
import firtez_dz as frz
import time

#start = time.time()
nx = 2048
ny = 2048
nz = 192
atm = frz.atm_model3D(nx,ny,nz)

T = np.fromfile("/dat/schmassmann/rempel/spot_12x8x12km_ng/temp.float", dtype = np.float32)
T = T.reshape(4096,256,4096).transpose(0,2,1) #reshape to (x,z,y) then transpose to (x,y,z)
atm.set_tem(T[:2048,2048:,:192]) # select x,y #from 0-2048 #
T = None
print('::temp saved',atm.tem.shape)

P = np.fromfile("/dat/schmassmann/rempel/spot_12x8x12km_ng/pres.float", dtype = np.float32)
P = P.reshape(4096,256,4096).transpose(0,2,1)
atm.set_pg(P[:2048,2048:,:192])
P = None
print('::pres saved')

m = np.fromfile("/dat/schmassmann/rempel/spot_12x8x12km_ng/magz.float", dtype = np.float32)
m = m.reshape(4096,256,4096).transpose(0,2,1)
atm.set_bx(m[:2048,2048:,:192])
m = None
m = np.fromfile("/dat/schmassmann/rempel/spot_12x8x12km_ng/magx.float", dtype = np.float32)
m = m.reshape(4096,256,4096).transpose(0,2,1)
atm.set_by(m[:2048,2048:,:192])
m = None
m = np.fromfile("/dat/schmassmann/rempel/spot_12x8x12km_ng/magy.float", dtype = np.float32)
m = m.reshape(4096,256,4096).transpose(0,2,1)
atm.set_bz(m[:2048,2048:,:192])
m = None
print('::mag saved')

v = np.fromfile("/dat/schmassmann/rempel/spot_12x8x12km_ng/velz.float", dtype = np.float32)
v = v.reshape(4096,256,4096).transpose(0,2,1)
atm.set_vx(v[:2048,2048:,:192])
v = None
v = np.fromfile("/dat/schmassmann/rempel/spot_12x8x12km_ng/velx.float", dtype = np.float32)
v = v.reshape(4096,256,4096).transpose(0,2,1)
atm.set_vy(v[:2048,2048:,:192])
v = None
v = np.fromfile("/dat/schmassmann/rempel/spot_12x8x12km_ng/vely.float", dtype = np.float32)
v = v.reshape(4096,256,4096).transpose(0,2,1)
atm.set_vz(v[:2048,2048:,:192])
v = None
print('::vel saved')

z = np.arange(192) * 8
z3d = np.zeros((nx,ny,nz))
z3d[:,:,:] = z[None,None,:]
atm.set_z(z3d)
print('::z saved')

atm.write_model('/dat/xenosh/atm_2048x2048x192_01.bin')
print('::atmos saved')

atm = None
#end = time.time()
#print(end - start)









