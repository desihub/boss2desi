import sys
import numpy as np
import fitsio
from numpy import linalg
from scipy.special import erf
import time
import scipy as sp
from scipy.interpolate import interp1d
from matplotlib import pyplot as pp
import iminuit
import traceback

def resolution(i,i0,sigma):
    ## This is the resolution in pixel units
    ## integrate the gaussian over 1 pixel

    res = (erf((i-i0+0.5)/sigma/sp.sqrt(2)) - erf((i-i0-0.5)/sigma/sp.sqrt(2)))*0.5*sigma
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

def newFitArc(arcfile,wave_new,arclines,fiber=None,debug=False,out=None,log=None):
    harc = fitsio.FITS(arcfile)
    flux = harc[0].read()
    nfib = flux.shape[0]
    wave = 10**harc[3].read()
    ok = sp.ones(nfib,dtype=bool)

    med = sp.median(flux.sum(axis=1))
    w = flux.sum(axis=1) < med/100.
    ok[w] = False

    wdisp = np.zeros([nfib,len(wave_new)])
    to = np.loadtxt(arclines,usecols=(0,))
    if fiber==None:
        fiber=range(nfib)
    else:
        fiber=[fiber]
    for fib in fiber:
        sys.stderr.write("fitting arc in fiber {}\n ".format(fib))
        index = np.arange(flux.shape[1])
        i = interp1d(wave[fib,:],index)
        w = (to>wave[fib,:].min()) & (to<wave[fib,:].max())
        try:
            wd,a,b = fitDisp(flux[fib,:],i(to[w]))
            wd = interp1d(wave[fib,:],wd,bounds_error=False,fill_value=wd.mean())
            wdisp[fib,:] = wd(wave_new)
            sys.stderr.wrote("mean(wdisp) in fib {} {}\n".format(fib,wdisp[fib,:].mean()))
            if debug:
                pp.figure(1)
                pp.plot(a,b)
                pp.plot(flux[fib,:])
                #for l in arclines:
                #    pp.plot(l+sp.zeros(2),flux[fib,:].max()*sp.arange(2),"k--")
                pp.show()
                pp.draw()
                x=raw_input()
                if x=="stop":
                    return wdisp,ok
                pp.clf()
        except:
            print "wdisp fit in fiber {} failed".format(fib)
            if debug:
                traceback.print_exc()
                return wdisp,ok
            wdisp[fib,:]=0.1
            ok[fib]=False
            if log is not None:
                log.write("failed wdisp fit in fiber {} \n".format(fib))

    if out is not None:
        fout=open(out,"w")
        for fib in range(wave.shape[0]):
            iok = 1
            if not ok[fib]:
                iok = 0
            fout.write("{} ".format(iok))
        fout.write("\n")
        for i in range(len(wave_new)):
            fout.write("{} ".format(wave_new[i]))
            for fib in range(wave.shape[0]):
                fout.write("{} ".format(wdisp[fib,i]))
            fout.write("\n")
        fout.close()
    return wdisp,ok
        
def fitDisp(flux,ilines,deg=2,tol=10):
    ''' 
    given a flux from the arc lamps fit sigmas
    '''
    index = sp.arange(len(flux))*1.
    p = sp.zeros(deg+1)
    imax = 1.*ilines.max()
    imin = 1.*ilines.min()

    dlam = abs(index[:,None]-ilines).min(axis=1)
    w=dlam<tol
    i0=index[w]
    f0=flux[w]

    def sigma(x,*p):
        s = x*0.
        u = (x-imin)/(imax-imin)
        for i in range(deg+1):
                s+=p[i]*u**i
        return sp.exp(s)

    def peaks(s):
        ## find the amplitudes of the peaks with a linear fit
        res = resolution(i0,ilines[:,None],s[:,None])
        D = (f0*res).sum(axis=1)
        M = res.dot(res.T)
        A = linalg.inv(M).dot(D)
        return A.T.dot(res)

    def chi2(*p):
        s = sigma(ilines,*p)
        res = peaks(s)
        ret = ((f0-res)**2).sum()
        return ret
    
    pinit={'a1': -0.341836969722409, 'a0': 0.5172259852497609, 'a2': 0.4141426841348048}
    pnames = ["a{}".format(i) for i in range(deg+1)]
    kwds = {p:pinit[p] for p in pnames}
    for p in pnames:
        kwds["error_"+p]=0.1
    kwds["limit_a0"]=(0,None)
    mig = iminuit.Minuit(chi2,forced_parameters=pnames,errordef=1,print_level=0,**kwds)
    mig.migrad()

    s0 = sigma(ilines,*[mig.values[p] for p in pnames])
    s1 = sigma(index,*[mig.values[p] for p in pnames])
    return s1,i0,peaks(s0)


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

def fitSkyLines(flux,ivar,ilines,tol=10):
    '''
    fit shift and dilation parameters to fix drifts in the wavelength solution from the arc
    '''
    
    nlines = len(ilines)
    nbins = len(flux)
    index = sp.arange(nbins)
    used_ilines = []
    fl = []
    iv = []
    x = []
    for i in range(nlines):
        w = abs(index-ilines[i]) < tol
        norm = ivar[w].sum()
        if norm==0:
            continue
        
        fl.append(flux[w])
        iv.append(ivar[w]/ivar[w].sum())
        x.append(index[w])
        used_ilines.append(ilines[i])
    fl = sp.array(fl)
    iv = sp.array(iv)
    x = sp.array(x)
    used_ilines = sp.array(used_ilines)
    nlines = len(used_ilines)
    ## require at least 4 lines to do the fits
    if len(used_ilines)<4:
        return index,0,0

    def peak(i,i0,sigma):
        f = sp.exp(-(i-i0)**2/2/sigma**2)
        xbar = (iv*f).sum(axis=1)
        x2bar = (iv*f**2).sum(axis=1)
        dbar = (iv*fl).sum(axis=1)
        dxbar= (iv*fl*f).sum(axis=1)
        M = [ sp.array([ [1,xbar[l]], [xbar[l],x2bar[l]] ]) for l in range(nlines)]
        C = [ linalg.inv(M[l]).dot([dbar[l],dxbar[l]]) for l in range(nlines)]
        C = sp.array(C)
        B = C[:,0,None]
        A = C[:,1,None]
        fit = B+A*f
        return fit

    def line_shift(i,epsilon,eta):
        return epsilon+(1+eta)*i

    def chi2(*p):
        p = sp.array(p)
        epsilon = p[0]
        eta = p[1]
        sigma = p[2:]
        peaks = sp.zeros((nlines,nbins))
        i = line_shift(x,epsilon,eta)
        peaks = peak(i,used_ilines[:,None],sigma[:,None])
        chisq = ((fl-peaks)**2*iv).sum()
        return ((fl-peaks)**2*iv).sum()

    pars=["epsilon","eta"]
    kwds={"epsilon":0.,"eta":0.,"error_epsilon":0.1,"error_eta":1e-4}
    kwds["limit_epsilon"]=(-2,2)
    kwds["limit_eta"]=(-1e-3,1e-3)
    for i in range(nlines):
        p = "sigma"+str(i)
        p0 = 1.
        dp0=0.5
        kwds[p]=p0
        kwds["error_"+p]=dp0
        kwds["limit_"+p]=(1./sp.sqrt(12),2.)
        pars.append(p)
#    kwds["fix_epsilon"]=True
#    kwds["fix_eta"]=True
    mig = iminuit.Minuit(chi2,forced_parameters=pars,errordef=1,print_level=0,**kwds)
    mig.migrad()
    epsilon = mig.values["epsilon"]
    eta = mig.values["eta"]
    sigma = [mig.values["sigma"+str(i)] for i in range(nlines)]
    sigma = sp.array(sigma)
    i = line_shift(index,epsilon,eta)
    return i,epsilon,eta

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
    ## round 4*std to get ndiag
    ndiag = int(4*std.max()+0.5)
    if ndiag%2==0:
        ndiag+=1

    print ndiag,std.max(),std.argmax()

    reso = sp.zeros([ndiag,nbins])
    for i in range(ndiag):
        offset = ndiag/2-i
        d = sp.diagonal(Q,offset=offset)
        if offset<0:
            reso[i,:len(d)] = d
        else:
            reso[i,nbins-len(d):nbins]=d

    t = time.time()
    print "spectro perfected in: ",t-t0
    return flux,ivar,reso

