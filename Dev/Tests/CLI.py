from functions import *

# INPUTS

imgdirpath = "Downsizing/guitar_origin_C001H001S0001_20220101_144452/"
refimgname = "guitar_origin_C001H001S0001_20220101_1444520001.tif"
angledata = "GuitarData.h5"

sqx = 0
sqy = 0

hnx = 1024
hny = 1024

patternx = 1.18533E-3
patterny = 1.18533E-3

F = 0.44

# END OF INPUTS

t0 = time.time()

imagelist  = sorted(list(filter(IS_IMAGE,os.listdir(imgdirpath))))
if (refimgname in imagelist):
    imagelist.remove(refimgname)

# imagelist=['guitar_origin_C001H001S0001_20220101_1444522950.tif']
imagelist = imagelist[2950-10:2950+10]

anglefile = h5py.File(angledata,'r')
angledata = anglefile['DATA']['{}'.format(sqy)]['{}'.format(sqx)]
D = angledata['D'][()]*1E-3
theta = angledata['Theta'][()]
phi = angledata['Phi'][()]

H = D*np.cos(phi)

# resultfile = h5py.File('({},{}).h5'.format(sqx,sqy),'w')
resultfile = h5py.File('test.h5','w')
udata = resultfile.create_dataset('u',shape=[len(imagelist),hny,hnx],dtype=np.double)
vdata = resultfile.create_dataset('v',shape=[len(imagelist),hny,hnx],dtype=np.double)
hdata = resultfile.create_dataset('h',shape=[len(imagelist),hny,hnx],dtype=np.double)
idata = resultfile.create_dataset('imagenames',shape=[len(imagelist)],dtype=h5py.string_dtype())

c1,c2,k1,k2,nx,ny = PROCESS_IMAGE(imgdirpath+refimgname)
scalex,scaley = GET_SCALE(H,F,theta,phi,k1,k2,patternx,patterny)
scalex = scalex*nx/hnx
scaley = scaley*nx/hny

resultfile.create_dataset('scalex',data=scalex)
resultfile.create_dataset('scaley',data=scaley)

A = INTGRAD2_A(hnx,hny,dx=scalex,dy=scaley,a11=False)

C0 = np.zeros([hny,hnx])

t1 = time.time()

print('\n[{}] Initialization complete'.format(TIMEC(t1-t0)))

for num,image in enumerate(imagelist):

    t2 = time.time()

    imagepath = os.path.join(imgdirpath,image)
    m1,m2,tp,tp,tp,tp = PROCESS_IMAGE(imagepath)
    u,v = GET_DISPLACEMENT(m1,m2,c1,c2,k1,k2)
    u = u[::(ny//hny),::(nx//hnx)]*scalex
    v = v[::(ny//hny),::(nx//hnx)]*scaley
    h,C0 = GET_HEIGHT(A,u,v,H,theta,scalex,scaley,C0=C0.flatten())

    tw = time.time()

    idata[num] = image
    udata[num] = u
    vdata[num] = v
    hdata[num] = h

    t3 = time.time()

    print(TIMEC(t3-tw))

    ela = t3-t0
    eta = t3-t1
    prg = (num+1)/len(imagelist)
    eta = eta*(1-prg)/prg
    prg = int(prg*100)
    print('[{}] Processed {}, Progress = {:3d}%, Elapsed {}, ETA {}'.format(TIMEC(t3-t2),image,prg,TIMEC(ela),TIMEC(eta)))

resultfile.close()

t4 = time.time()
print('[{}] Calculation completed\n'.format(TIMEC(t4-t0)))