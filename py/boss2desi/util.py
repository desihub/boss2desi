import numpy as np
import fitsio
from numpy import linalg
from scipy.special import erf
import time
import scipy as sp
from scipy.interpolate import interp1d
from matplotlib import pyplot as pp

def resolution(i,i0,sigma):
    ## This is the resolution in pixel units
    ## integrate the gaussian over 1 pixel
    res = (erf((i-i0+0.5)/sigma) - erf((i-i0-0.5)/sigma))*0.5*sigma
    
    return res

def fitPeaks(index,flux,npeaks,i0,w0,debug=False):

    tol = 5 ##  width in pixels around maximum where to fit the peak

    if npeaks == 0:
        return

    imax = flux.argmax()
    index_max = index[imax]
    w=(index>index_max-tol) & (index<index_max+tol)
    i1 = np.average(index[w],weights = flux[w])

    wdisp = 4.*np.arange(1000)/1000.+1/sp.sqrt(12)
    re = resolution(index[w],i1,wdisp[:,None])
    norm = re.sum(axis=1)
    wn = norm>0
    re[wn,:]/=norm[wn,None]
    chi2 = ((flux[w]/flux[w].sum()-re)**2).sum(axis=1)
    iwdisp = chi2.argmin()

    if debug:
        pp.figure(1)
        pp.plot(index[w],flux[w]/flux[w].sum())
        pp.plot(index[w],re[iwdisp])
        pp.show()
        pp.draw()
        x=raw_input()
        pp.clf()

    i0.append(i1)
    w0.append(wdisp[iwdisp])

    ## now mask and findPeak again
    fitPeaks(index[~w],flux[~w],npeaks-1,i0,w0,debug=debug)

def fitArc(arcfile,camera,lambda_out=None):
    h = fitsio.FITS(arcfile)

    flux = h[0].read()
    lam = 10**h[3].read()
    nfibers = flux.shape[0]
    nbins = flux.shape[1]
    if lambda_out is not None:
        nbins = len(lambda_out)
    
    wdisp = np.zeros([nfibers,nbins])
    nfib = flux.shape[0]

    nlines = 25
    if camera[0]=="r":
        nlines = 15
    deg = 2
    for fib in range(nfib):
        l0 = []
        w0 = []

        l = lam[fib,:]
        f = flux[fib,:]
        mask = ((l<3640) | (l > 3670)) & ((l<4330) | (l>4351))
        i = sp.arange(len(l))
        i = i[mask]
        l = l[mask]
        f = f[mask]
        fitPeaks(i,f,nlines,l0,w0)
        l0 = np.array(l0)
        w0 = np.array(w0)

        M = np.zeros([deg+1,nlines])
        for a in range(deg+1):
            M[a,:]=l0**a

        md = (M*w0).sum(axis=1)
        A = np.zeros([deg+1,deg+1])
        for a in range(deg+1):
            for b in range(a,deg+1):
                A[a,b] = (M[a,:]*M[b,:]).sum()
                A[b,a] = A[a,b]
        
        A = linalg.inv(A)
        p = A.dot(md)
        wd = l*0
        for a in range(deg+1):
            wd+=p[a]*i**a
        
        if lambda_out is not None:
            wd = interp1d(l,wd,bounds_error=False,fill_value=0)
            wdisp[fib,:]=wd(lambda_out)
        else:
            wdisp[fib,mask]=wd
        
    return lam,wdisp,l0,w0

def spectro_perf(fl,iv,re):
    t0 = time.time()
    ## compute R and F
    R = sp.sqrt(iv)*re
    F = R.dot(sp.sqrt(iv)*fl)
    R = R.dot(R.T)

    ## diagonalize
    d,s = linalg.eigh(R)

    ## invert
    d_inv = d*0
    w = d>d.max()*1e-10
    d[~w]=0
    d_inv[w]=1/d[w]
    IR = s.dot(d_inv[:,None]*s.T)
    F = IR.dot(F)
    Q = s.dot(sp.sqrt(d[:,None])*s.T)

    norm = Q.sum(axis=1)
    w=norm>0
    Q[w,:] = Q[w,:]/norm[w,None]
    flux = Q.dot(F)
    ivar = norm**2

    ## find the variance of the Q matrices
    ## ndiag ~ 3*sqrt(var)

    nbins = len(flux)
    index = sp.arange(nbins)
    mean = (Q*index).sum(axis=1)
    ## add at least the variance of a pixel in case Q is a delta function
    var = (Q*(index-mean[:,None])**2).sum(axis=1)+1./12
    std = sp.sqrt(var)
    ## round 3*std to get ndiag
    ndiag = int(3*std.max()+0.5)
    if ndiag%2==0:
        ndiag+=1
    
    print ndiag,std.max(),std.argmax()

    reso = sp.zeros([ndiag,nbins])
    for i in range(ndiag):
        reso[i,:nbins-abs(i-ndiag/2)] = sp.diagonal(Q,offset=i-ndiag/2)

    t = time.time()
    print "spectro perfected in: ",t-t0
    return flux,ivar,reso

