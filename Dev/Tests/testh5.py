import numpy as np
import h5py
import matplotlib.pyplot as pl

f = h5py.File('test.h5','r')
h = f['h'][:]
u = f['u'][:]
v = f['v'][:]
sx = f['scalex'][()]
sy = f['scaley'][()]

print(sx,sy)
exit()

g = h5py.File("GuitarData.h5",'r')['DATA']['0']['0']
H = g['H'][()]
alpha = g['Alpha'][()]
theta = g['Theta'][()]

T = 0.5*np.sin(2*theta)
H = H*np.cos(alpha)

C = np.gradient(h,sy,axis=1) + (T*h/H) + (0.5*v/H)
# C = C*H/T

# C = np.gradient(h,sx,axis=2) + (0.5*u/H)

for i in range(h.shape[0]):
    pl.figure()
    pl.title('Frame = {}'.format(i))
    im = pl.pcolormesh(C[i],cmap='jet',vmin=np.amin(C[:]),vmax=np.amax(C[:]))
    pl.colorbar(im)
    pl.show()

# D = np.mean(C,axis=(-2,-1))
# E = np.ones_like(C)
# for i in range(h.shape[0]):
#     E[i] = (C[i]-D[i])**2
# E = np.mean(E,axis=(-2,-1))
# E = E**0.5

# pl.figure()
# pl.plot(D)
# pl.plot(E)
# pl.show()
