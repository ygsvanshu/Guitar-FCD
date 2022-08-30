from functions import *
from matplotlib import rc
from matplotlib.animation import FFMpegWriter

t0 = time.time()

rc('font',**{'family':'STIXGeneral','serif':['stix']})
rc('text', usetex=True)
rc('font',**{'size':13})

# INPUTS

sqx = 0
sqy = 0

# END OF INPUTS

F = h5py.File('({},{}).h5'.format(sqx,sqy),'r')

tsteps = F['h'].shape[0]
ndt = int(np.log10(tsteps))+1
mdt = min(ndt,5)+1

ny = F['h'].shape[1]
nx = F['h'].shape[2]
scalex = F['scalex'][()]
scaley = F['scaley'][()]

X = np.linspace(0,(nx-1)*scalex,nx)
Y = np.linspace(0,(ny-1)*scaley,ny)
X, Y = np.meshgrid(X, Y)

hmin = np.amin(F['h'][:])*1E5
hmax = np.amax(F['h'][:])*1E5

fig, ax = pl.subplots(subplot_kw={"projection": "3d"},constrained_layout=True)
ax.set_zlim(hmin,hmax)
ax.set_xlim(0,(nx-1)*scalex)
ax.set_ylim(0,(ny-1)*scaley)
ax.zaxis.set_rotate_label(False)
ax.set_xlabel('x [m]')
ax.set_ylabel('y [m]')
ax.set_zlabel('h [m]'+r'$\times 10^{-5}$',rotation=90)
ax.view_init(30,225)

surf = ax.plot_surface(X, Y, F['h'][0]*1E5, cmap='jet', vmin=hmin, vmax=hmax, linewidth=0, antialiased=False)
pl.colorbar(surf)

writer = FFMpegWriter(fps=120)

with writer.saving(fig,'({},{}).mp4'.format(sqx,sqy),dpi=200):
    t1 = time.time()
    for i in range(tsteps):
        t2 = time.time()
        surf.remove()
        surf = ax.plot_surface(X, Y, F['h'][i]*1E5, cmap='jet', vmin=hmin, vmax=hmax, linewidth=0, antialiased=False)
        ax.set_title('Frame {:0{w1}d}/{}, time = {:{w2}.4f} sec'.format(i+1,tsteps,(i+1)/10000,w1=ndt,w2=mdt))
        pl.draw()
        writer.grab_frame()
        prg = (i+1)/tsteps
        t3 = time.time()
        ela = t3-t0
        eta = t3-t1
        eta = eta*(1-prg)/prg
        prg = int(prg*100)
        print('[{}] Processed {:{w}d}, Progress = {:3d}%, Elapsed {}, ETA {}'.format(TIMEC(t3-t2),i+1,prg,TIMEC(ela),TIMEC(eta),w=ndt))

t4 = time.time()
print('[{}] Movie generated\n'.format(TIMEC(t4-t0)))