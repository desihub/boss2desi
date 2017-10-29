from __future__ import print_function
import numpy as np
from scipy.interpolate import interp1d
from .util import svd_spectro_perf

from multiprocessing import Pool,cpu_count

class spectrum:
    def __init__(self,wave,flux,ivar,mask,wdisp):
        self.flux = flux
        self.wave = wave
        self.ivar = ivar
        self.mask = mask
        self.wdisp = wdisp

class frame:
    def __init__(self, plate, mjd, camera, expid, flavor, spectra):
        self.camera = camera
        self.expid = expid
        self.flavor = flavor
        self.plate = plate
        self.mjd = mjd

        self.spectra = np.array(spectra)

def _wrap_pool(spec):
    wave_new = spec.wave_new
    nbins = len(spec.wave)
    i = np.arange(nbins)
    lam_to_i = interp1d(spec.wave,i)
    i0 = lam_to_i(wave_new)
    di2 = (i-i0[:,None])**2
    re = np.exp(-di2/2/spec.wdisp**2)
    re /= re.sum(axis=1)[:,None]

    ## select wave_new pixels "close" to an unmasked pixel
    ## an unmasked pixel:
    w = spec.ivar >0
    if w.sum()==0:
        print("DEBUG: skipping spectrum with no data")
        return wave_new*0,wave_new*0,wave_new*0+1,np.zeros((spec.ndiag,len(wave_new)))
    mask = abs((i0-i[w,None])/spec.wdisp[w,None]).min(axis=0) < 0.5

    try:
        flux,ivar,R = svd_spectro_perf(spec.flux,spec.ivar,re)
    except:
        print("DEBUG: svd didn't converge, skipping")
        return wave_new*0,wave_new*0,wave_new*0+2,np.zeros((spec.ndiag,len(wave_new)))
    ivar *= mask

    ndiag = spec.ndiag
    nbins = R.shape[1]
    reso = np.zeros([ndiag,nbins])
    for i in range(ndiag):
        offset = ndiag/2-i
        d = np.diagonal(R,offset=offset)
        if offset<0:
            reso[i,:len(d)] = d
        else:
            reso[i,nbins-len(d):nbins]=d

    return flux,ivar,~mask,reso

def resample(fra,wave_new,ndiag,ncpu=None):

    ## pack wave_new and ndiag into spectra to pass to the pool of multi-threads
    for s in fra.spectra:
        s.wave_new = wave_new
        s.ndiag = ndiag    

    if ncpu is None:
        ncpu = cpu_count()

    if ncpu > 1:
        pool = Pool(ncpu)
        print("INFO: Starting resampling on {} processors".format(ncpu))

        res =  pool.map(_wrap_pool,fra.spectra)
        pool.close()
    elif ncpu == 1:
        res =  map(_wrap_pool,fra.spectra)
        
    spectra = [spectrum(s.wave_new,fl,iv,ma,R) for s,(fl,iv,ma,R) in zip(fra.spectra,res)]

    #plate, mjd, camera, expid, flavor, spectra
    return frame(fra.plate,fra.mjd,fra.camera,fra.expid,fra.flavor,spectra)

