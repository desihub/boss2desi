import scipy as sp
import fitsio
import boss2desi.fibermap

class brick:
    def __init__(self,fl,iv,ma,wave,wave_new,wdisp,camera,fibermap,sky=None,wave_wdisp=None):

        self.camera = camera

        la = wave
        la_wdisp = wave_wdisp
        wd = wdisp

        ## ndiag = 3 wdisp
        self.ndiag = 2*int(3*wdisp.max())+1
        print "ndiag=",self.ndiag,wdisp.max()

        ## turn off bits coming from sky subtraction:
        ma = ma & (~2**22) & (~2**23)

        ## turn off bits coming from flat-fielding
        ma = ma & (~2**17)

        if not sky is None:
            fl+=sky
        nspec = fl.shape[0]
        nbins = len(wave_new)
        lam = wave_new
        flux = sp.zeros([nspec,nbins])
        ivar = sp.zeros([nspec,nbins])
        mask = sp.zeros([nspec,nbins],dtype=sp.uint32)
        re = sp.zeros([nspec,self.ndiag,nbins])

        self.fm = None
        if fibermap is not None:
            desi_fibermap = boss2desi.fibermap.boss_fibermap(fibermap,self.camera)
            self.fm = desi_fibermap.fm
            self.fm_names = desi_fibermap.fm_names


        for i in range(nspec):
            j = sp.searchsorted(la[i,:],lam)
            w=j>=len(la[i,:])
            j[w]-=1
            flux[i,:] = (la[i,j]-lam)*fl[i,j-1]*iv[i,j-1] + (lam-la[i,j-1])*fl[i,j]*iv[i,j]
            mask[i,:] = ma[i,j-1] & ma[i,j]
            norm = (la[i,j]-lam)*iv[i,j-1] + (lam-la[i,j-1])*iv[i,j]
            ivar[i,:] = norm**2

            norm_ivar = (iv[i,j-1]*(la[i,j]-lam)**2 + iv[i,j]*(lam-la[i,j-1])**2)
            w=norm_ivar>0
            ivar[i,w]/=norm_ivar[w]

            w=(iv[i,j-1]==0) | (iv[i,j]==0)
            norm[w]=0
            flux[i,w]=0
            ivar[i,w]=0

            w=norm>0
            flux[i,w]/=norm[w]

            wdisp=(la_wdisp[i,j]-lam)*wd[i,j-1] + (lam-la_wdisp[i,j-1])*wd[i,j]
            wdisp/=la_wdisp[i,j]-la_wdisp[i,j-1]
            re[i,:,:] = sp.exp(-(sp.arange(self.ndiag)-self.ndiag/2)[:,None]**2/2./wdisp**2)
            re[i,:,:]/=sp.sum(re[i,:,:],axis=0)

        self.lam = lam
        self.flux = flux
        self.ivar = ivar
        self.mask = mask
        self.re = re
    
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



