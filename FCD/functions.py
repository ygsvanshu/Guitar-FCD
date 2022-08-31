def PRINT_ERROR(string):
    print('\033[91m'+'[ERROR !!!!] '+string+'\033[0m')
    return()

def PRINT_SUCCESS(string):
    print('\033[92m'+'[SUCCESS !!] '+string+'\033[0m')
    return()

def PRINT_WARNING(string):
    print('\033[93m'+'[WARNING !!] '+string+'\033[0m')
    return()

def ASK_YN(question):
    question += ' (Y/N): '
    answer = 'A'
    while(answer not in ['Y','y','N','n']):
        answer = input(question)
    if ((answer=='Y')or(answer=='y')): return(True)
    if ((answer=='N')or(answer=='n')): return(False)

def IMPORT(module,alias=None):

    module0 = module.split('.')[0]
    module1 = module.split('.')[-1]

    try:
        if alias: 
            globals()[alias] = importlib.import_module(name=module)
        else: 
            globals()[module1] = importlib.import_module(name=module)
    except:
        PRINT_WARNING('Cannot import the module: {}'.format(module0))
        install = ASK_YN('Try to install module: {}'.format(module0))
        if install:
            os.system('pip3 install --upgrade {}'.format(module0))
            os.system('conda update {}'.format(module0))

    try:
        if alias: 
            globals()[alias] = importlib.import_module(name=module)
        else: 
            globals()[module] = importlib.import_module(name=module)
    except:
        PRINT_ERROR('Cannot import the module: {}'.format(module0))
        exit()

    return()

##[Vanshu] System modules (default in all python3 installations)
from logging import raiseExceptions
import os
import sys
import time
import importlib

##[Vanshu] Additional modules (need to be checked for the current python3 installation)
modules = ['numpy','h5py','scipy.io','scipy.ndimage','scipy.fft','scipy.sparse.linalg','scipy.sparse','PIL.Image','matplotlib.pyplot','matplotlib.colors']
aliases = ['np',None,'spio','sndi','sfft','sslg','sspr',None,'pl','cl']

##[Vanshu] Import all additional modules in the list
for mod,ali in zip(modules,aliases):
    IMPORT(module=mod,alias=ali)

def TIMEC(s):

    h = int(s/3600)
    s = s-(h*3600)
    m = int(s/60)
    s = s-(m*60)
    # s = int(s)

    return('{:02d}:{:02d}:{:06.3f}'.format(h,m,s))

def TIMED(s):

    h = int(s/3600)
    s = s-(h*3600)
    m = int(s/60)
    s = s-(m*60)
    s = int(s)

    return('{:02d}:{:02d}:{:02d}'.format(h,m,s))

def ANGLE1(idx,kx,ky):
    x = kx[idx[1]]
    y = ky[idx[0]]
    a = np.arctan2(y,x)
    return(a)

def IS_IMAGE(filename):
    extension = filename.split('.')[-1]
    if (extension in ['bmp','dds','eps','ico','jpg','jpeg','png','tga','tif','tiff']):
        return(True)
    else:
        return(False)

##[Vanshu] Function extracts carrier frequencies and carriers from image
def PROCESS_IMAGE(img_name,gau_s=3,uni_s=9,pks_n=4,exp_s=2,plots=False):

    ##[Vanshu] These inputs need to be integers
    pks_n = int(pks_n)
    exp_s = int(abs(exp_s))

    ##[Vanshu] Opening the image
    img = Image.open(img_name)
    # print('Image format = ',img.format, ', Image size = ',img.size,', Image mode = ',img.mode)
    img = np.array(img)
    igx = img.shape[1]
    igy = img.shape[0]

    ##[Vanshu] Performing the FFT and getting the frequencies
    fft = sfft.fftshift(sfft.fft2(img))
    kx  = 2*np.pi*sfft.fftshift(sfft.fftfreq(img.shape[1]))
    ky  = 2*np.pi*sfft.fftshift(sfft.fftfreq(img.shape[0]))

    frl = np.real(fft)
    fim = np.imag(fft)

    ##[Vanshu] Isolating the highest 'pks_n' number of peaks
    fgf = np.fabs(fim)
    fgf = sndi.gaussian_filter(fgf,sigma=gau_s)             ##[Vanshu] Gaussian filter to smooth it out
    fgf[fgf<sndi.maximum_filter(fgf,size=uni_s)] = 0        ##[Vanshu] Finding the local maxima
    fgf[ky==0,kx==0] = 0                                    ##[Vanshu] Ignoring the peak at kx=0 and ky=0 
    fgf = np.argpartition(fgf,-pks_n,axis=None)[-pks_n:]    ##[Vanshu] Taking indices of the highest 'pks_n' number of peaks in flattened form

    ##[Vanshu] Converting them from flattened indices to 2D indices and sorting them in clockwise manner
    pks = []
    for pk in fgf:
        pks += [np.unravel_index(pk,fim.shape)]
    pks = sorted(pks,key=lambda a:ANGLE1(a,kx,ky))
    pks = np.array(pks)
    
    ##[Vanshu] Finding the average radius of circular filter around the peaks
    dst = 0
    for i,pk in enumerate(pks):
        pk2 = pks[i]
        pk1 = pks[i-1]
        dkx = kx[pk2[1]]-kx[pk1[1]]
        dky = ky[pk2[0]]-ky[pk1[0]]
        dst = dst + ((dkx**2)+(dky**2))**0.5
    dst = dst/pks.shape[0]
    rad = dst/2

    ##[Vanshu] Choosing k1 and k2 based on maximum cross product
    imx = 0
    jmx = 0
    crs = 0
    for i in range(pks_n):
        for j in range(pks_n):
            kx1 = kx[pks[i][1]]
            kx2 = kx[pks[j][1]]
            ky1 = ky[pks[i][0]]
            ky2 = ky[pks[j][0]]
            km1 = ((kx1**2) + (ky1**2))**0.5
            km2 = ((kx2**2) + (ky2**2))**0.5
            krs = abs(((kx1*ky2) - (kx2*ky1))/(km1*km2))
            if (krs>crs):
                crs = krs
                imx = i
                jmx = j

    ##[Vanshu] For the first carrier 
    msk = np.zeros_like(fim)
    ppx = kx[pks[imx][1]]
    ppy = ky[pks[imx][0]]
    kkx = (kx-ppx)**2
    kky = (ky-ppy)**2
    kkx,kky = np.meshgrid(kkx,kky)
    dst = (kkx + kky)**0.5
    if (exp_s>0):
        msk = np.cos(0.5*np.pi*((dst/rad)**exp_s))**2   ##[Vanshu] Fancy cosine type fall-off (use exp_s>0)
        msk[dst>rad]=0
    else:
        msk[dst<=rad]=1                                 ##[Vanshu] Simple binary cutoff (use exp_s=0)

    k1 = [kx[pks[2][1]],ky[pks[2][0]]]                  ##[Vanshu] First carrier wavenumber
    c1 = sfft.ifft2(sfft.ifftshift(fft*msk))            ##[Vanshu] First carrier

    ##[Vanshu] For the second carrier 
    msk = np.zeros_like(fim)
    ppx = kx[pks[jmx][1]]
    ppy = ky[pks[jmx][0]]
    kkx = (kx-ppx)**2
    kky = (ky-ppy)**2
    kkx,kky = np.meshgrid(kkx,kky)
    dst = (kkx + kky)**0.5
    if (exp_s>0):
        msk = np.cos(0.5*np.pi*((dst/rad)**exp_s))**2   ##[Vanshu] Fancy cosine type fall-off (use exp_s>0)
        msk[dst>rad]=0
    else:
        msk[dst<=rad]=1                                 ##[Vanshu] Simple binary cutoff (use exp_s=0)

    k2 = [kx[pks[3][1]],ky[pks[3][0]]]                  ##[Vanshu] Second carrier wavenumber
    c2 = sfft.ifft2(sfft.ifftshift(fft*msk))            ##[Vanshu] Second carrier

    ##[Vanshu] Plots for debugging purposes. Turning this on will cause a loss of performance
    if plots:

        fig,axs = pl.subplots(1,3,figsize=(14,4),constrained_layout='True')
        
        ax = axs[0]
        ax.set_title('Image')
        axs[0].imshow(img,cmap='binary')
        ax.set_aspect(1)
        ax.set_xlabel(r'$x$')
        ax.set_ylabel(r'$y$')
        
        ax = axs[1]
        ax.set_title('Real')
        im = ax.pcolor(kx,ky,np.fabs(frl),norm=cl.SymLogNorm(linthresh=1e-3),cmap='jet')
        pl.colorbar(im,ax=ax)
        ax.set_aspect(1)
        ax.set_xlabel(r'$k_x$')
        ax.set_ylabel(r'$k_y$')

        ax = axs[2]
        ax.set_title('Imaginary')
        im = ax.pcolor(kx,ky,np.fabs(fim),norm=cl.SymLogNorm(linthresh=1e-3),cmap='jet')
        pl.colorbar(im,ax=ax)
        for i,pk in enumerate(pks):
            aa = np.linspace(0,2*np.pi,37)
            xx = kx[pk[1]] + rad*np.cos(aa)
            yy = ky[pk[0]] + rad*np.sin(aa)
            ax.plot(kx[pk[1]],ky[pk[0]],linestyle='',marker='+',color='k')
            ax.plot(xx,yy,color='k',linestyle='-')
            ax.text(kx[pk[1]],ky[pk[0]],i+1)
        ax.set_aspect(1)
        ax.set_xlabel(r'$k_x$')
        ax.set_ylabel(r'$k_y$')

        pl.savefig(img_name.split('.')[0]+'.png',dpi=300,bbox_inches='tight')
        pl.close(fig)

    return(c1,c2,k1,k2,igx,igy)
  
##[Vanshu] Function computes Laplacian
def LAP(f):

    M = f.shape[1]
    N = f.shape[0]

    p = 2*np.pi*sfft.fftshift(sfft.fftfreq(M))
    q = 2*np.pi*sfft.fftshift(sfft.fftfreq(N))

    p,q = np.meshgrid(p,q)

    g = sfft.fftshift(sfft.fft2(f))
    g = g*((p**2)+(q**2))
    g = sfft.ifft2(sfft.ifftshift(g))
    g = -np.real(g)

    return(g)

##[Vanshu] Function computes inverse Laplacian
def ILAP(f):

    M = f.shape[1]
    N = f.shape[0]

    p = 2*np.pi*sfft.fftshift(sfft.fftfreq(M))
    q = 2*np.pi*sfft.fftshift(sfft.fftfreq(N))

    p[p==0]=np.inf
    q[q==0]=np.inf

    p,q = np.meshgrid(p,q)

    g = sfft.fftshift(sfft.fft2(f))
    g = g/((p**2)+(q**2))
    g = sfft.ifft2(sfft.ifftshift(g))
    g = -np.real(g)

    return(g)

##[Vanshu] Correction for wrapped phase based on 
##[Vanshu] Volkov, Vyacheslav V., and Yimei Zhu. "Deterministic phase unwrapping in the presence of noise." Optics letters 28.22 (2003): 2156-2158.
##[Vanshu] Jeught SVd, Sijbers J, Dirckx JJJ. Fast Fourier-Based Phase Unwrapping on the Graphics Processing Unit in Real-Time Imaging Applications. Journal of Imaging. 2015; 1(1):31-44. https://doi.org/10.3390/jimaging1010031
def PHASE_WRAP(f):

    k = np.cos(f)*LAP(np.sin(f)) - np.sin(f)*LAP(np.cos(f)) - LAP(f)
    k = np.round((1/(2*np.pi))*ILAP(k))
    f = f + (2*np.pi*k)

    return(f)

##[Vanshu] Function to get displacements (solving 2x2 matrix per pixel)
def GET_DISPLACEMENT(m1,m2,c1,c2,k1,k2):

    kx1  = k1[0]
    ky1  = k1[1]

    kx2  = k2[0]
    ky2  = k2[1]

    phi1 = np.imag(np.log(c1*np.conj(m1)))  ##[Vanshu] Computing the phase modulation for first carrier
    phi2 = np.imag(np.log(c2*np.conj(m2)))  ##[Vanshu] Computing the phase modulation for second carrier

    phi1 = PHASE_WRAP(phi1)                 ##[Vanshu] Correcting the wrapped phase for first carrier
    phi2 = PHASE_WRAP(phi2)                 ##[Vanshu] Correcting the wrapped phase for second carrier

    det  = (kx1*ky2) - (kx2*ky1)            ##[Vanshu] Computing the determinant for inversion of 2x2 matrix

    u    = ((ky2*phi1) - (ky1*phi2))/det    ##[Vanshu] Computing the horizontal displacement
    v    = ((kx1*phi2) - (kx2*phi1))/det    ##[Vanshu] Computing the vertical displacement

    return(u,v)

def INTGRAD2_A(nx,ny,dx=1,dy=1):

    if ((not isinstance(nx,int))or(not isinstance(ny,int))):
        raise ValueError('nx and ny must be integers')

    if (min(nx,ny)<2):
        raise ValueError('nx and ny must be at least 2')

    npdx = np.array(dx)
    if ((npdx.size!=1)and(npdx.shape!=nx)):
        raise ValueError('dx is not a scalar or of length nx')

    npdy = np.array(dy)
    if ((npdy.size!=1)and(npdy.shape!=ny)):
        raise ValueError('dy is not a scalar or of length ny')

    if (npdx.size==1):
        dxp = dx*np.ones([nx])
        dxm = dx*np.ones([nx])
    else:
        dxp = np.pad((dx[1:]-dx[:-1]),(0,1),'edge')
        dxm = np.pad((dx[1:]-dx[:-1]),(1,0),'edge')

    if (npdy.size==1):
        dyp = dy*np.ones([ny])
        dym = dy*np.ones([ny])
    else:
        dyp = np.pad((dy[1:]-dy[:-1]),(0,1),'edge')
        dym = np.pad((dy[1:]-dy[:-1]),(1,0),'edge')

    k = 0

    if (nx%2==0):
        if (ny%2==0):
            I = np.zeros([6*nx*ny+4])
            J = np.zeros([6*nx*ny+4])
            C = np.zeros([6*nx*ny+4])

            I[0] = 0
            J[0] = ((ny//2)-1)*nx + (nx//2) - 1
            C[0] = 0.25

            I[1] = 0
            J[1] = ((ny//2)-1)*nx + (nx//2)
            C[1] = 0.25

            I[2] = 0
            J[2] = (ny//2)*nx + (nx//2) - 1
            C[2] = 0.25

            I[3] = 0
            J[3] = (ny//2)*nx + (nx//2)
            C[3] = 0.25

            k = 4
        else:
            I = np.zeros([6*nx*ny+2])
            J = np.zeros([6*nx*ny+2])
            C = np.zeros([6*nx*ny+2])

            I[0] = 0
            J[0] = (ny//2)*nx + (nx//2) - 1
            C[0] = 0.5

            I[1] = 0
            J[1] = (ny//2)*nx + (nx//2)
            C[1] = 0.5

            k = 2
    else:
        if (ny%2==0):
            I = np.zeros([6*nx*ny+2])
            J = np.zeros([6*nx*ny+2])
            C = np.zeros([6*nx*ny+2])

            I[0] = 0
            J[0] = ((ny//2)-1)*nx + (nx//2)
            C[0] = 0.5

            I[1] = 0
            J[1] = (ny//2)*nx + (nx//2)
            C[1] = 0.5

            k = 2
        else:
            I = np.zeros([6*nx*ny+1])
            J = np.zeros([6*nx*ny+1])
            C = np.zeros([6*nx*ny+1])

            I[0] = 0
            J[0] = (ny//2)*nx + (nx//2)
            C[0] = 1.0

            k = 1

    l = 1
    for j in range(ny):
        for i in range(nx):
            
            I[k  ] = l
            I[k+1] = l
            I[k+2] = l

            ic = j*nx + i
            im = (ic - 1)%(nx*ny)
            ip = (ic + 1)%(nx*ny)

            J[k  ] = im
            J[k+1] = ic
            J[k+2] = ip

            d_p = dxp[i]
            d_m = dxm[i]
            
            if (i==0):
                C[k  ] = 0.0
                C[k+1] = -1/d_p
                C[k+2] = 1/d_p
            elif(i==nx-1):
                C[k  ] = -1/d_m
                C[k+1] = 1/d_m
                C[k+2] = 0.0
            else:
                C[k  ] = -d_p/(d_m*(d_m+d_p))
                C[k+1] = (d_p-d_m)/(d_m*d_p)
                C[k+2] = d_m/(d_p*(d_m+d_p))

            k = k+3
            l = l+1
    
    for j in range(ny):
        for i in range(nx):

            I[k  ] = l
            I[k+1] = l
            I[k+2] = l

            jc = j*nx + i
            jm = (jc - nx)%(nx*ny)
            jp = (jc + nx)%(nx*ny)

            J[k  ] = jm
            J[k+1] = jc
            J[k+2] = jp

            d_p = dyp[j]
            d_m = dym[j]
            
            if (j==0):
                C[k  ] = 0.0
                C[k+1] = -1/d_p
                C[k+2] = 1/d_p
            elif(j==ny-1):
                C[k  ] = -1/d_m
                C[k+1] = 1/d_m
                C[k+2] = 0.0
            else:
                C[k  ] = -d_p/(d_m*(d_m+d_p))
                C[k+1] = (d_p-d_m)/(d_m*d_p)
                C[k+2] = d_m/(d_p*(d_m+d_p))

            k = k+3
            l = l+1

    A = sspr.csc_matrix((C,(I,J)),shape=(2*nx*ny+1,nx*ny))

    AT = (A.transpose()).tocsc()
    AA = (AT.dot(A)).tocsc()

    SS = sspr.linalg.factorized(AA)

    return(SS,AT)

def INTGRAD2_B(SS,AT,fx,fy,h00=0):

    if ((fx.ndim!=2)or(fy.ndim!=2)):
        raise ValueError('fx and fy must be 2d arrays')

    if (fx.shape!=fy.shape):
        raise ValueError('fx and fy must be the same sizes.')

    if (AT.shape[0]!=fx.size):
        raise ValueError('Matrix A not compatible for the given fx and fy; re-generate matrix A')

    if ((not isinstance(h00,(int,float)))or(not np.isfinite(h00))or(not np.isreal(h00))):
        raise ValueError('h00 must be a finite scalar numeric variable.') 

    B = np.concatenate((np.array([h00]),fx.flatten(),fy.flatten()))
    B = AT.dot(B)
    C = SS(B)
    C = C.reshape(fx.shape)

    return(C)

def GET_SCALE(H,F,theta,alpha,k1,k2,pattern_x,pattern_y):
    
    kx1 = k1[0]
    ky1 = k1[1]

    kx2 = k2[0]
    ky2 = k2[1]

    if (abs(ky1)<abs(ky2)):
        Gxx = (2*np.pi)/(((kx1**2)+(ky1**2))**0.5)
        Gyy = (2*np.pi)/(((kx2**2)+(ky2**2))**0.5)
    else:
        Gxx = (2*np.pi)/(((kx2**2)+(ky2**2))**0.5)
        Gyy = (2*np.pi)/(((kx1**2)+(ky1**2))**0.5)

    factor = 1 - (H*np.cos(alpha))/(F*np.cos(theta))

    scalex = (pattern_x/Gxx)*(factor/np.cos(alpha))
    scaley = (pattern_y/Gyy)*(factor)

    return(scalex,scaley)

def GET_HEIGHT(SS,AT,u,v,H,theta,scalex=1,scaley=1,h00=0):

    nx = u.shape[1]
    ny = u.shape[0]

    T = 0.5*np.sin(2*theta)

    x  = np.array([scalex*np.arange(nx)])
    x  = x - 0.5*(x[0,-1])
    y  = np.transpose(np.array([scaley*np.arange(ny)]))
    y  = y - 0.5*(y[-1,0])

    ig = -np.exp(x*T/H)/(2*H)

    fu = ig*u
    fv = ig*v

    fh = INTGRAD2_B(SS,AT,fu,fv,h00=h00)
    h  = fh/np.exp(x*T/H)

    return(h)