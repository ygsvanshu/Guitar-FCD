import numpy as np
import h5py

d1 = np.load('NPYs/day1.npy')
d2 = np.load('NPYs/day2.npy')
d3 = np.load('NPYs/day3.npy')
d4 = np.load('NPYs/day4.npy')
d5 = np.load('NPYs/day5.npy')

readme = open('README.txt').read()

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

ff = h5py.File('GuitarData.h5','w')
rr = ff.create_dataset('README',data=readme)
gg = ff.create_group('DATA')

for i,r in enumerate(rows):
    
    cc = columns[i]
    rw = gg.create_group('{}'.format(r))

    if (r==-4):
        day = 5
        data = d5[:12]
    elif (r==-3):
        day = 5
        data = d5[12:]
    elif(r==-2):
        day = 1
        data = d1
    elif(r==-1):
        day = 2
        data = d2
    elif(r==0):
        day = 2
        data = d2
    elif(r==1):
        day = 3
        idx = [0] + [*range(1,d3.shape[0],2)]
        data = d3[idx]
    elif(r==2):
        day = 3
        idx = [*range(0,d3.shape[0],2)][1:]
        data = d3[idx]
    elif(r==3):
        day = 4
        idx = [0] + [*range(1,d4.shape[0],2)]
        data = d4[idx]
    elif(r==4):
        day = 4
        idx = [*range(0,d4.shape[0],2)][1:]
        data = d4[idx]
    else:
        print('SHITE')

    print('R=',r,data.shape)

    for j,c in enumerate(cc):

        dd = data[j]
        print('   C=',c,dd.shape)
        co = rw.create_group('{}'.format(c))
        co.create_dataset('Day',data=day)
        co.create_dataset('LightPoint',data=[dd[0],dd[1]])
        co.create_dataset('LightLine',data=[[dd[2],dd[3]],[dd[4],dd[5]]])
        co.create_dataset('CameraLine',data=[[dd[6],dd[7]],[dd[8],dd[9]]])
        co.create_dataset('Phi',data=dd[10])
        co.create_dataset('D',data=dd[11])
        co.create_dataset('Theta',data=dd[12])

ff.close()    