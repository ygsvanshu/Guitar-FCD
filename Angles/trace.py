import os
import numpy as np
import h5py
import matplotlib.pyplot as pl

# trace = np.genfromtxt('Trace/Trace.csv',delimiter=',')
# centroid = np.mean(trace,axis=0)
# angle = np.arctan2((trace[:,1]-centroid[1]),(trace[:,0]-centroid[0]))
# trace = trace[np.argsort(angle)]
# np.savetxt('Trace/Trace.dat',trace)

rows = [*range(-4,5,1)]
columns = [
    [*range(-7,5,1)],
    [*range(-8,5,1)],
    [*range(-8,4,1)],
    [*range(-8,6,1)],
    [*range(-8,6,1)],
    [*range(-8,5,1)],
    [*range(-7,5,1)],
    [*range(-2,5,1)],
    [*range(-1,5,1)]
]

trace = np.genfromtxt('Trace/Trace.dat',skip_header=3)
plane = np.genfromtxt('Trace/Trace.dat',max_rows=2)
centroid = np.mean(trace,axis=0)
trace = trace-centroid
plane = plane-centroid

angle = np.arctan2(abs(plane[0,0]-plane[1,0]),abs(plane[0,1]-plane[1,1]))

trace2 = np.zeros_like(trace)
plane2 = np.zeros_like(plane)

trace3 = np.zeros_like(trace)
plane3 = np.zeros_like(plane)

trace2[:,0] = -trace[:,1]
trace2[:,1] = -trace[:,0]

plane2[:,0] = -plane[:,1]
plane2[:,1] = -plane[:,0]

trace3[:,0] = trace2[:,0]*np.cos(angle) - trace2[:,1]*np.sin(angle)
trace3[:,1] = trace2[:,0]*np.sin(angle) + trace2[:,1]*np.cos(angle)

plane3[:,0] = plane2[:,0]*np.cos(angle) - plane2[:,1]*np.sin(angle)
plane3[:,1] = plane2[:,0]*np.sin(angle) + plane2[:,1]*np.cos(angle)

centroid3 = np.mean(trace3,axis=0)

hline = np.mean(plane3[:,1])
rline = np.amax(trace3[:,0])

origin = [rline-210,hline+180]

fig,axs = pl.subplots(figsize=(10,8))
axs.plot(trace3[:,0],trace3[:,1],'k-')
axs.plot(plane3[:,0],plane3[:,1],'g-')
# axs.plot(centroid3[0],centroid3[1],'m*')
for i in range(16):
    axs.axvline(rline-30*i,linewidth=0.5,color='gray')
for i in range(13):
    axs.axhline(hline+30*i,linewidth=0.5,color='gray')
for i,r in enumerate(rows):
    for j,c in enumerate(columns[i]):
        axs.text(origin[0]+30*c+15,origin[1]+30*r+10,'({},{})'.format(c,r),horizontalalignment='center')
axs.set_aspect(1)
axs.set_xlabel('mm')
axs.set_ylabel('mm')
pl.savefig('Trace/Grid.png')

