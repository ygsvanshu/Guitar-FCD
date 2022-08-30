from functions import *
import multiprocessing as mp
from functools import partial

def CALCULATE(num,imglist,imgdir,c1,c2,k1,k2,A,H,theta,scalex,scaley):

    t0 = time.time()

    image = imglist[num]
    imagepath = os.path.join(imgdir,image)
    m1,m2,tp,tp,tp,tp = PROCESS_IMAGE(imagepath)
    u,v = GET_DISPLACEMENT(m1,m2,c1,c2,k1,k2)
    h = GET_HEIGHT(A,u,v,H,theta,scalex,scaley,h11=0,a11=False)

    t1 = time.time()
    return(u,v,h,image,t1-t0)

def RUN(imgdir,H,F,theta,phi,patternx,patterny):

    t0 = time.time()

    imglist = sorted(filter(IS_IMAGE,os.listdir(imgdir)))
    refimg = imglist.pop(0)
    
    c1,c2,k1,k2,nx,ny = PROCESS_IMAGE(os.path.join(imgdir,refimg))
    scalex,scaley = GET_SCALE(H,F,theta,phi,k1,k2,patternx,patterny)
    A = INTGRAD2_A(nx,ny,dx=scalex,dy=scaley,a11=False)

    resultfile = h5py.File('test.h5','w')
    resultfile.create_dataset('u',shape=[len(imglist),nx,ny],dtype=np.double)
    resultfile.create_dataset('v',shape=[len(imglist),nx,ny],dtype=np.double)
    resultfile.create_dataset('h',shape=[len(imglist),nx,ny],dtype=np.double)
    resultfile.close()

    CALC = partial(CALCULATE,imglist=imglist,imgdir=imgdir,c1=c1,c2=c2,k1=k1,k2=k2,A=A,H=H,theta=theta,scalex=scalex,scaley=scaley)

    print('Initialised')

    with mp.Pool() as pool:
        resultfile = h5py.File('test.h5','r+')
        results = pool.imap(CALC,range(len(imglist)))
        for num,result in enumerate(results):
            resultfile['u'][num]=result[0]
            resultfile['v'][num]=result[1]
            resultfile['h'][num]=result[2]
            print("Done for {} in {}".format(result[-2],TIMEC(result[-1])))
        resultfile.close()
    pool.close()
    pool.join()

    t1 = time.time()

    print('\n[{}] Computation complete'.format(TIMEC(t1-t0)))

if (__name__=='__main__'):
    imgdir = sys.argv[1] #'/Users/vanshu/Documents/Acoustic/plus1minus3_C001H001S0001_20211215_164624-10'
    sqx = sys.argv[2]
    sqy = sys.argv[3] 
    patternx = 1.18533E-3
    patterny = 1.18533E-3

    anglefile = h5py.File('GuitarData.h5','r')
    angledata = anglefile['DATA']['{}'.format(sqy)]['{}'.format(sqx)]
    D = angledata['D'][()]*1E-3
    theta = angledata['Theta'][()]
    phi = angledata['Phi'][()]
    F = 0.44
    H = D*np.cos(phi)
    RUN(imgdir,H,F,theta,phi,patternx,patterny)