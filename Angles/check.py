import os
import numpy as np
import h5py
from PIL import Image
import matplotlib.pyplot as pl


days = [1,2,3,4,5]
data = h5py.File('GuitarData.h5','r')

if not os.path.exists('Verification'): os.mkdir('Verification')

sx = 297#*(2**0.5)
sy = 420#*(2**0.5)

mm = 0.1/2.54

for day in days:

    im = Image.open("PNGs/guitarday{}.png".format(day))

    if (day==5): 
        sx = 297*(2**0.5)
        sy = 420*(2**0.5)

    fig, ax = pl.subplots(figsize=(sx*mm,sy*mm))
    ax.imshow(im,cmap='binary_r',extent=[0,sx,0,sy])

    for r in list(data['DATA'].keys()):
        for c in list (data['DATA'][r].keys()):
            dt = data['DATA'][r][c]
            if (dt['Day'][()]==day):
                l1 = ax.plot(dt['LightLine'][:,0],dt['LightLine'][:,1],color='r',marker='',linestyle='--',linewidth=0.5)  # LightLine
                l2 = ax.plot(dt['CameraLine'][:,0],dt['CameraLine'][:,1],color='g',marker='',linestyle='--',linewidth=0.5)  # CameraLine
                p0 = ax.plot(dt['LightPoint'][0],dt['LightPoint'][1],color='b',marker='+',linestyle='',linewidth=0.5)  # LightPoint
    pl.axis('off')
    pl.margins(0,0)
    pl.savefig('Verification/day{}.png'.format(day),dpi=300,bbox_inches='tight',pad_inches=0)
    pl.close(fig)