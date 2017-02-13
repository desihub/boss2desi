import scipy as sp
import fitsio



class brick:
    ndiag = 11
    def __init__(self,fl,iv,ma,wave,wave_new,wdisp,camera,fibermap,sky=None,wave_wdisp=None):

        self.camera = camera

        iv *= (ma==0)
        la = wave
        la_wdisp = wave_wdisp
        wd = wdisp
        if not sky is None:
            fl+=sky
        nspec = fl.shape[0]
        nbins = len(wave_new)
        lam = wave_new
        flux = sp.zeros([nspec,nbins])
        ivar = sp.zeros([nspec,nbins])
        mask = sp.zeros([nspec,nbins],dtype=sp.uint32)
        re = sp.zeros([nspec,self.ndiag,nbins])

        self.read_fibermap(fibermap)
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
    
    def read_fibermap(self,fibermap):
        ## get fibermap names
        fm_names = fibermap.get_colnames()
        nspec = len(fibermap[fm_names[0]][:])
        ## create a fiber map and add all columns from boss' fibermap
        fm = []
        for i in fm_names:
            fm.append(fibermap[i][:])

        ## add fiber number
        fm_names.append("FIBER")
        fib_num = sp.arange(nspec,dtype=int)+1
        if self.camera == 'b2' or self.camera == 'r2':
            fib_num += 500
        fm.append(fib_num)

        ## add filter names
        fm_names.append("FILTER")
        fm.append(sp.array([["SDSS_U","SDSS_G","SDSS_R","SDSS_I","SDSS_Z"]]*nspec))

        ## correct magnitudes
        for i,obj in enumerate(fm[fm_names.index('MAG')]):
            ## avoid NA and SKY
            ## sky fibers don't have magnitudes
            if fm[fm_names.index('OBJTYPE')][i].strip() != 'SKY'\
                    and fm[fm_names.index('OBJTYPE')][i].strip() != 'NA':
                ## fix only positive fluxes...
                w = fm[fm_names.index('CALIBFLUX')][i] > 0
                fm[fm_names.index("MAG")][i][w] = \
                        22.5-2.5*sp.log10(fm[fm_names.index('CALIBFLUX')][i][w])

        ## change names of standard stars to desi convention
        for i,obj in enumerate(fm[fm_names.index("OBJTYPE")]):
            fm[fm_names.index("OBJTYPE")][i] = obj.replace("SPECTROPHOTO_STD","STD")

        self.fm_names = fm_names
        self.fm = fm

    def export(self,fout="brick.fits.gz"):
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



