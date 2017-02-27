import scipy as sp
from numpy import linalg
from scipy import sparse
from scipy.sparse.linalg import lobpcg
from scipy.interpolate import interp1d

import fitsio
import boss2desi.fibermap
from boss2desi import util

class brick:
    def __init__(self,fl,iv,ma,wave,wave_new,wdisp,camera,fibermap):

        self.camera = camera


        ## turn off bits coming from sky subtraction:
        ma = ma & (~2**22) & (~2**23)

        ## turn off bits coming from flat-fielding
        ma = ma & (~2**17)

        nspec = fl.shape[0]
        nbins = len(wave_new)
        flux = sp.zeros([nspec,nbins])
        ivar = sp.zeros([nspec,nbins])
        mask = sp.zeros([nspec,nbins],dtype=sp.uint32)
        re = []

        self.fm = None
        if fibermap is not None:
            desi_fibermap = boss2desi.fibermap.boss_fibermap(fibermap,self.camera)
            self.fm = desi_fibermap.fm
            self.fm_names = desi_fibermap.fm_names

        index = sp.arange(len(wave[0,:]))
        for fib in range(nspec):
            print fib
            i_wave = interp1d(wave[fib,:],index)
            wlam = (wave[fib,:]>wave_new.min()) & (wave[fib]<wave_new.max())
            res = util.resolution(index[wlam],i_wave(wave_new[:,None]),wdisp[fib,:,None])
            norm = res.sum(axis=1)
            res/=norm[:,None]

            ## check that the centers of the resolution are within
            ## two pixels of a good pixel
            w = iv[fib,:]>0
            if w.sum()==0:
                re.append(sp.zeros([2,nbins]))
                continue
            centers = res.dot(index[wlam])
            wgood = abs(centers-index[wlam & w,None]).min(axis=0)<2
            if wgood.sum()==0:
                re.append(sp.zeros([2,nbins]))
                continue
            f,i,r = util.spectro_perf(fl[fib,wlam],iv[fib,wlam],res[wgood,:])
            flux[fib,wgood]=f
            ivar[fib,wgood]=i
            reso = sp.zeros([r.shape[0],nbins])
            reso[:,wgood]=r
            re.append(reso)
            #f,i,r = util.spectro_perf(fl[fib,wlam],iv[fib,wlam],res)
            #flux[fib,:]=f
            #ivar[fib,:]=i
            #re.append(r)

        ndiags = sp.array([r.shape[0] for r in re])
        w = ivar.sum(axis=1)>0
        ndiag = max(ndiags[w])
        print "max in diag {} in fiber {}".format(ndiag,sp.argmax(ndiags))
        self.re = sp.zeros([nspec,ndiag,nbins])
        for fib in range(len(re)):
            nd = re[fib].shape[0]
            self.re[fib,(ndiag-nd)/2:(ndiag+nd)/2]=re[fib]

        self.lam = wave_new
        self.flux = flux
        self.ivar = ivar
        self.mask = mask
    
    def export(self,fout):
        hlist=[]
        hlist.append({"name":"CAMERA","value":self.camera,"comment":" Spectograph Camera"})

        brick = fitsio.FITS(fout,"rw",clobber=True)
        brick.write(self.flux,extname="FLUX",header=hlist)
        brick.write(self.ivar,extname="IVAR")
        brick.write(self.mask,extname="MASK")
        brick.write(self.lam,extname="WAVELENGTH")
        brick.write(self.re,extname="RESOLUTION")
        brick.write(self.fm,names=self.fm_names,extname="FIBERMAP")
        brick.close()

        '''
        hlist = []
        hlist.append({"name":"CRVAL1","value":7089.0,"comment":" Starting wavelength [Angstroms]"})
        hlist.append({"name":"CDELT1","value":0.2,"comment":" Wavelength step [Angstroms]"})
        hlist.append({"name":"AIRORVAC","value":"vac","comment":" Vacuum wavelengths"})
        hlist.append({"name":"LOGLAM","value":0,"comment":" linear wavelength steps, not log10"})
        hlist.append({"name":"SIMFILE","value":"alpha-3/20150107/simspec-00000003.fits","comment":" Input simulation file"})
        hlist.append({"name":"VSPECTER","value":"0.0.0","comment":" TODO: Specter version"})
        hlist.append({"name":"EXPTIME","value":1000.0,"comment":" Exposure time [sec]"})
        hlist.append({"name":"RDNOISE","value":2.9,"comment":" Read noise [electrons]"})
        hlist.append({"name":"FLAVOR","value":"science","comment":" Exposure type (arc, flat, science)"})
        hlist.append({"name":"SPECMIN","value":0,"comment":" First spectrum"})
        hlist.append({"name":"SPECMAX","value":1,"comment":" Last spectrum"})
        hlist.append({"name":"NSPEC","value":2,"comment":" Number of spectra"})
        hlist.append({"name":"WAVEMIN","value":7089.0,"comment":" First wavelength [Angstroms]"})
        hlist.append({"name":"WAVEMAX","value":7110.5,"comment":" Last wavelength [Angstroms]"})
        hlist.append({"name":"WAVESTEP","value":0.5,"comment":" Wavelength step size [Angstroms]"})
        hlist.append({"name":"SPECTER","value":"0.1.dev1","comment":" https://github.com/sbailey/specter"})
        hlist.append({"name":"IN_PSF","value":"psf-r-bidon.fits","comment":" Input spectral PSF"})
        hlist.append({"name":"IN_IMG","value":"pix-r0-bidon.fits","comment":" Input image"})
        '''



