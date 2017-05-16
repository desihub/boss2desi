import sys
import numpy as np
import fitsio
from scipy import linalg
from scipy.special import erf
import time
import scipy as sp
from scipy.interpolate import interp1d
from matplotlib import pyplot as pp
import iminuit
import traceback

def resolution(i,i0,sigma,dpix=sp.sqrt(2)):
    ## This is the resolution in pixel units
    ## integrate the gaussian over 1 pixel
    res = (erf((i-i0+dpix)/sigma/sp.sqrt(2)) - erf((i-i0-dpix)/sigma/sp.sqrt(2)))*0.5*sigma
    ##res = (erf((i-i0+.5)/sigma/sp.sqrt(2)) - erf((i-i0-.5)/sigma/sp.sqrt(2)))*0.5*sigma
    ##res = sp.exp(-(i-i0)**2/2/sigma**2)
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

def newFitArc(arcfile,wave_new,arclines,fiber=None,debug=False,out=None,log=None,deg=None,deg_bb=None,tol=None):
    harc = fitsio.FITS(arcfile)
    flux = harc[0].read()
    ivar = harc[1].read()
    nfib = flux.shape[0]
    wave = 10**harc[3].read()
    ok = sp.ones(nfib,dtype=bool)

    med = sp.median(flux.sum(axis=1))
    w = flux.sum(axis=1) < med/100.
    ok[w] = False

    wdisp = np.zeros([nfib,len(wave_new)])
    dpix  = np.zeros(nfib)
    to = np.loadtxt(arclines,usecols=(0,))
    index = np.arange(flux.shape[1])
    if fiber==None:
        fiber=range(nfib)
    for fib in fiber:
        sys.stderr.write("fitting arc in fiber {}\n ".format(fib))
        i = interp1d(wave[fib,:],index)
        w = (to>wave_new.min()) & (to<wave_new.max())
        try: 
            t0 = time.time()
            sys.stderr.write("fitDisp\n")
            wd,a,b,pars,chi2,ndf,dpix[fib] = fitDisp(flux[fib,:],ivar[fib,:],i(to[w]),deg=deg,log=log,deg_bb=deg_bb,tol=tol)
            sys.stderr.write("fit Disp done in {}".format(time.time()-t0))
            print "fit Disp in ",time.time()-t0

            wd = interp1d(wave[fib,:],wd,bounds_error=False,fill_value=wd.mean())
            wdisp[fib,:] = wd(wave_new)
            if log is not None:
                log.write("wdisp({}) ".format(fib))
                for p in pars:
                    log.write("{} ".format(p))
                log.write("\n")
                log.flush()
            sys.stderr.write("mean(wdisp) in fib {} {}, dpix {}\n".format(fib,wdisp[fib,:].mean(),dpix[fib]))
            if debug:
                sys.stderr.write("chi2: {}, ndf: {} ".format(chi2,ndf))
                pp.figure(1)
                pp.plot(wave[fib,a],b)
                pp.plot(wave[fib,:],flux[fib,:]*(ivar[fib,:]>0))
                pp.grid()
                for l in to:
                    pp.plot(l+sp.zeros(2),flux[fib,:].max()*sp.arange(2),"k--")
                pp.show()
                pp.draw()
                x=raw_input()
                if x=="stop":
                    return wdisp,ok
                pp.clf()
        except:
            sys.stdout.write("wdisp fit in fiber {} failed\n".format(fib))
            if debug:
                traceback.print_exc()
                return wdisp,ok,dpix
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
    return wdisp,ok,dpix
        
def fitDisp(flux,ivar,ilines,tol=10,deg=2,log=None,p0=None,deg_bb=3):
    ''' 
    given a flux from the arc lamps fit sigmas
    '''
    if deg is None:
        deg=2
    if deg_bb is None:
        deg_bb=3
    if tol == None:
        tol = 10
    nbins = len(flux)
    index = sp.arange(nbins)
    imax = 1.*nbins
    imin = 0.
    nlines=len(ilines)

    ## nodes at the zeros of the chebyshev polynomials to minimize ringing?
    ##node = sp.cos(sp.pi*(2*sp.arange(deg+1,dtype=float)+1)/2/(deg+1))
    node = -1+2*sp.arange(deg+1,dtype=float)/deg
    pol = sp.ones([deg+1,len(ilines)])
    u = (2.*ilines-(imin+imax))/(imax-imin)
    for i in range(deg+1):
        for j in range(deg+1):
            if j==i:continue
            pol[i]*=(u-node[j])/(node[i]-node[j])

    dlam = abs(index[:,None]-ilines).min(axis=1)
    w=dlam<tol
    i0=index[w]
    f0=flux[w]
    eta=0.0
    iv=ivar[w]/((eta*f0)**2*ivar[w]+1)

    def sigma(p):
        s=p.dot(pol)
        return sp.exp(s)

    def peaks(s,fl,dpix=sp.sqrt(2)):
        ## find the amplitudes of the peaks with a linear fit
        res = resolution(i0,ilines[:,None],s[:,None],dpix=dpix)
        ivres = sp.sqrt(iv)*res
        D = (sp.sqrt(iv)*fl*ivres).sum(axis=1)
        M = ivres.dot(ivres.T)
        A = linalg.inv(M).dot(D)

        return A.T.dot(res)

    def broadband(p):
        ret = p[0]+i0*0
        for i in range(1,deg_bb+1):
            ret += p[i]*i0**i
        return ret

    def chi2(*p):
        s = sigma(sp.array(p[:deg+1]))
        bb = broadband(p[deg+1:])
        res = peaks(s,f0-bb,dpix=p[-1])+bb
        ret = (iv*(f0-bb-res)**2).sum()
        return ret
    pinit={}
    if p0 is None:
        p0 = sp.ones(deg+1)
    for i in range(deg+1):
        pinit['a'+str(i)]=p0[i]
    pnames = ["a{}".format(i) for i in range(deg+1)]
    kwds = {p:pinit[p] for p in pnames}
    for p in pnames:
        kwds["error_"+p]=0.1
        kwds["limit_"+p]=(-0.5*sp.log(12.),sp.log(3))

    for i in range(deg_bb+1):
        p="bb_{}".format(i)
        pnames.append(p)
        kwds[p]=0.
        kwds["error_"+p]=0.1

    pnames.append("dpix")
    kwds["dpix"]=0.5
    kwds["error_dpix"]=0.1*sp.sqrt(2)
    kwds["limit_dpix"]=(0.,3)

    mig = iminuit.Minuit(chi2,forced_parameters=pnames,errordef=1,print_level=0,**kwds)
    t0 = time.time()
    mig.migrad()
    sys.stderr.write("mig in {}\n".format(time.time()-t0))
    dpix=mig.values["dpix"]

    pvals = sp.array([mig.values[p] for p in pnames])
    s0 = sigma(pvals[:deg+1])
    pol = sp.ones([deg+1,nbins])
    u = (2.*index-(imin+imax))/(imax-imin)
    for i in range(deg+1):
        for j in range(deg+1):
            if j==i:continue
            pol[i]*=(u-node[j])/(node[i]-node[j])
    s1 = sigma(pvals[:deg+1])
    bb = broadband(pvals[deg+1:])
    return s1,i0,peaks(s0,f0-bb,dpix)+bb,mig.values.values(),mig.fval,(iv>0).sum()-len(pnames),dpix

def fitArcLines(flux,ilines,tol=10,log=None,p0=None,ivar=None):
    ''' 
    given a flux from the arc lamps fit sigmas
    '''
    nbins = len(flux)
    index = sp.arange(nbins)
    nlines=len(ilines)

    dlam = abs(index[:,None]-ilines).min(axis=1)
    w=dlam<tol
    i0=index[w]
    f0=flux[w]
    if ivar==None:
        iv = f0*0+1
    else:
        iv = ivar[w]

    def peaks(s):
        ## find the amplitudes of the peaks with a linear fit
        res = resolution(i0,ilines[:,None],s[:,None])
        D = (f0*res).sum(axis=1)
        M = res.dot(res.T)
        A = linalg.inv(M).dot(D)
        return A.T.dot(res)

    def sigma(*p):
        return sp.array(p)

    def chi2(*p):
        s = sigma(p)
        res = peaks(s)
        ret = (iv*(f0-res)**2).sum()
        print sp.mean(s),ret
        return ret

    pinit={}
    if p0 is None:
        p0 = sp.ones(nlines)
    pnames = ["sigma{}".format(i) for i in range(nlines)]
    for i in range(nlines):
        pinit['sigma'+str(i)]=p0[i]
    kwds = {p:pinit[p] for p in pnames}
    for p in pnames:
        kwds["error_"+p]=0.1
        kwds["limit_"+p]=(1/sp.sqrt(12),3)
    mig = iminuit.Minuit(chi2,forced_parameters=pnames,errordef=1,print_level=3,**kwds)
    mig.migrad()

    s0 = sp.array([mig.values[p] for p in pnames])
    return i0,peaks(s0),mig.values.values()


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
        try:
            C = [ linalg.inv(M[l]).dot([dbar[l],dxbar[l]]) for l in range(nlines)]
        except:
            stop
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

def fitSkyLinesGlobally(flux,ivar,ilines,tol=10,deg_epsilon=2,deg_eta=1,log=None):
    '''
    fit shift and dilation parameters to fix drifts in the wavelength solution from the arc
    assume shift, dilatation and sigma are quadratic with fiber number
    '''
    
    fibers = sp.arange(flux.shape[0])
    index = sp.arange(flux.shape[1])
    ## works for 500 fibers and degree 2
    nodes_eps = fibers[0]+(fibers[-1]-fibers[0])*sp.arange(deg_epsilon+1,dtype=float)/(deg_epsilon-1)
    pol_ep = sp.ones([deg_epsilon+1,len(fibers)])
    for i in range(deg_epsilon+1):
        for j in range(deg_epsilon+1):
            if j!=i:
                pol_ep[i]*=(fibers-nodes_eps[j])/(nodes_eps[i]-nodes_eps[j])

    def epsilon(p):
        return (sp.array(p)[:,None]*pol_ep).sum(axis=0)

    ## works for 500 fibers and degree 1
    def eta(fib,p):
        et=p[0]*(499.-fib)/499.
        et+=p[1]*fib/499.

        return et
    
    ## fit constant and amplitude for a given peak
    def peak(x,i0,sigma,fl,iv):
        f = sp.exp(-(x-i0)**2/2/sigma**2)
        xbar = (iv*f).sum()
        x2bar = (iv*f**2).sum()
        dbar = (iv*fl).sum()
        dxbar= (iv*fl*f).sum()
        M = sp.zeros([2,2])
        C = sp.zeros(2)
        M[0,0]=iv.sum()
        M[1,0]=xbar
        M[0,1]=xbar
        M[1,1]=x2bar
        C[0]=dbar
        C[1]=dxbar

        try:
            C = linalg.inv(M).dot(C)
        except:
            return fl*0
        A = C[0]
        B = C[1]
        fit = A+B*f
        return fit

    def chi2(*p):
        p_ep = p[0:nep]
        p_et = p[nep:nep+net]
        sigma = p[nep+net]
        ep = epsilon(p_ep)
        et = eta(fibers,p_et)
        chi=0
        for fib in fibers:
            il = ilines[fib]
            wall = abs(index-il[:,None])<tol
            for i,w in enumerate(wall):
                fit=peak(ep[fib]+(1+et[fib])*index[w],il[i],sigma,flux[fib,w],ivar[fib,w])
                chi += ((flux[fib,w]-fit)**2*ivar[fib,w]).sum()
        return chi

    nep=deg_epsilon+1
    pars=[]
    for i in range(nep):
        pars.append("a{}_ep".format(i))
    net=deg_eta+1
    for i in range(net):
        pars.append("a{}_et".format(i))
    pars.append("sigma")
    kwds={}
    for p in pars:
        kwds[p]=0.
        kwds["error_"+p]=0.01
        kwds["limit_"+p]=(-2,2)

    for i in range(net):
        kwds["fix_a{}_et".format(i)]=True

    kwds["sigma"]=1.
    kwds["error_sigma"]=0.1
    kwds["limit_sigma"]=(1./sp.sqrt(12),3)
    mig = iminuit.Minuit(chi2,forced_parameters=pars,errordef=1,print_level=3,**kwds)
    mig.migrad()

    pvals = [mig.values[p] for p in pars]
    ep = epsilon(pvals[0:nep])
    et = eta(fibers,pvals[nep:nep+net])
    sigma=pvals[nep+net]

    peaks = sp.zeros(flux.shape)
    index = ep[:,None]+(1+et[:,None])*index
    for fib in fibers:
        il = ilines[fib]
        wall = abs(index[fib]-il[:,None])<tol
        for i,w in enumerate(wall):
            i0=ilines[fib,i]
            peaks[fib,w]=peak(index[fib,w],i0,sigma,flux[fib,w],ivar[fib,w])

    return index,ep,et,peaks

def spectro_perf(fl,iv,re,tol=1e-3,log=None,ndiag_max=27):
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
    #w=d>1e-20
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

    ## ndiag is such that the sum of the diagonals > 1-tol
    ## and at most ndiag_max

    ndiag=1
    for i in range(Q.shape[0]):
        imin = i-ndiag/2
        if imin<0:imin=0
        imax = i+ndiag/2
        if imax>=Q.shape[1]:
            imax = Q.shape[1]

        frac=Q[i,imin:imax].sum()
        while frac<1-tol:
            ndiag+=2
            imin = i-ndiag/2
            if imin<0:imin=0
            imax = i+ndiag/2
            if imax>=Q.shape[1]:imax = Q.shape[1]
            frac = Q[i,imin:imax].sum()

    if ndiag>ndiag_max:
        log.write("WARNING, reducing ndiag {} to {}".format(ndiag,ndiag_max))
        ndiag=ndiag_max
    nbins = Q.shape[1]
    reso = sp.zeros([ndiag,nbins])
    for i in range(ndiag):
        offset = ndiag/2-i
        d = sp.diagonal(Q,offset=offset)
        if offset<0:
            reso[i,:len(d)] = d
        else:
            reso[i,nbins-len(d):nbins]=d

    t = time.time()
    sys.stdout.write("spectro perfected in: {} \n".format(t-t0))
    if log is not None:
        log.write("\n final ndiag: {}\n".format(ndiag))
        log.write("spectro perfected in: {} \n".format(t-t0))
        log.flush()
    return flux,ivar,reso


def svd_spectro_perf(fl,iv,re,log=None):
    t0 = time.time()
    ## compute R and F
    R = sp.sqrt(iv)*re
    R = R.T
    F = sp.sqrt(iv)*fl

    ## svd decomposition
    u,s,vt = linalg.svd(R)
    one = linalg.diagsvd(s*0+1,R.shape[0],R.shape[1])
    s = linalg.diagsvd(s,R.shape[0],R.shape[1])

    flux = vt.T.dot(one.T.dot(u.T.dot(F)))
    Q = vt.T.dot(sp.sqrt(s.T.dot(s)).dot(vt))

    norm = Q.sum(axis=1)
    w=norm>0
    Q[w,:] = Q[w,:]/norm[w,None]
    flux[w]/=norm[w] 
    ivar = norm**2

    t = time.time()
    sys.stdout.write("spectro perfected in: {} \n".format(t-t0))
    if log is not None:
        log.write("spectro perfected in: {} \n".format(t-t0))
    return flux,ivar,Q
'''
    ## ndiag is such that the sum of the diagonals > 1-tol
    ndiag=1
    for i in range(Q.shape[0]):
        imin = i-ndiag/2
        if imin<0:imin=0
        imax = i+ndiag/2
        if imax>=Q.shape[1]:
            imax = Q.shape[1]

        frac=Q[i,imin:imax].sum()
        while frac<1-tol:
            ndiag+=2
            imin = i-ndiag/2
            if imin<0:imin=0
            imax = i+ndiag/2
            if imax>=Q.shape[1]:imax = Q.shape[1]
            frac = Q[i,imin:imax].sum()

    if ndiag>ndiag_max:
        if log is not None:
            log.write("WARNING, reducing ndiag {} to {}".format(ndiag,ndiag_max))
        ndiag=ndiag_max
    nbins = Q.shape[1]
    reso = sp.zeros([ndiag,nbins])
    for i in range(ndiag):
        offset = ndiag/2-i
        d = sp.diagonal(Q,offset=offset)
        if offset<0:
            reso[i,:len(d)] = d
        else:
            reso[i,nbins-len(d):nbins]=d
'''


def convert_air_to_vacuum(air_wave) :
    ## copied over from specex
    sigma2 = (1e4/air_wave)**2
    fact = 1. +  5.792105e-2/(238.0185 - sigma2) +  1.67917e-3/( 57.362 - sigma2)
    vacuum_wave = air_wave*fact

    return vacuum_wave
