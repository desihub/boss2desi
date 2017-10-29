from __future__ import print_function
import fitsio
import numpy as np
from numpy.polynomial.legendre import Legendre
from numpy.polynomial.chebyshev import Chebyshev
from desispec.resolution import Resolution
from os.path import dirname
from os import makedirs

from .frame import spectrum,frame


def parse_coeff(x,coeff,row):

    func = row["FUNC"][0]
    if func == "legendre":
        func = Legendre
    elif func == "chebyshev":
        func = Chebyshev

    xmin = row["XMIN"][0]
    xmax = row["XMAX"][0]
    try:
        xjumplo = row["XJUMPLO"][0]
        xjumphi = row["XJUMPHI"][0]
        xjumpval = row["XJUMPVAL"][0]

        t = (x-xjumplo)/(xjumphi-xjumplo)
        jump = (t<0) + t*(0<t<1) + (t>1)
        x = x + jump*xjumpval
    except:
        pass

    xmid = (xmax+xmin)/2
    mu = 2*(x-xmid)/(xmax-xmin)

    return func(coeff)(mu)

def get_frames(spplatename,fibers=None,flatdir=None):

    spplate = fitsio.FITS(spplatename)
    path = dirname(spplatename)
    head = spplate[0].read_header()
    file_exps = []
    file_flat = []
    flat_exp = []
    for nexpc in ["B1","B2","R1","R2"]:
        cardname = "NEXP_{}".format(nexpc)
        if cardname in head:
            nexp = head[cardname]
            for i in range(nexp):
                iexp = len(file_exps)+1
                str_iexp = str(iexp)
                if iexp < 10:
                    str_iexp = '0'+str_iexp
                c = head["EXPID{}".format(str_iexp)][0:2]
                file_exps.append(path+"/spFrame-"+head["EXPID{}".format(str_iexp)][:11]+".fits.gz")
                file_flat.append(path+"/spFlat-"+c+head["EXPID{}".format(str_iexp)][11:20]+".fits.gz")
                flat_exp.append(c+head["EXPID{}".format(str_iexp)][11:20])
                

    spplate.close()
    flat_exp = np.unique(flat_exp)

    frames = []
    if flatdir is not None:
        plate = head["PLATEID"]
        mjd = head["MJD"]
        flavor='flat'
        for fe in flat_exp:
            flatfile = flatdir+"/"+str(plate)+"/flat-{}.fits".format(fe)
            print("INFO: reading flat file: {}".format(flatfile))
            h=fitsio.FITS(flatfile)
            flux = h[0].read()
            mask = h[2].read()
            ivar = h[1].read()
            loglam = h[3].read()
            wdisp = h[4].read()
            if fibers is not None:
                flux=flux[fibers]
                ivar=ivar[fibers]
                mask=mask[fibers]
                loglam=loglam[fibers]
                wdisp=wdisp[fibers]
            spectra = []
            for i in range(flux.shape[0]):
                spectra.append(spectrum(10**loglam[i],flux[i],ivar[i],mask[i],wdisp[i]))
                
            frames.append(frame(plate,mjd,fe[:2],int(fe[3:]),flavor,spectra))


    for i,fexp in enumerate(file_exps):
        fra = get_frame(fexp,file_flat[i],fibers=fibers)
        frames.append(fra)

    return frames

def get_frame(fexp,flat,fibers=None):
    print("INFO: reading frame {}".format(fexp))
        
    h = fitsio.FITS(fexp)
    spFlat = fitsio.FITS(flat)
    c = h[0].read_header()["CAMERAS"].strip()
    expid = h[0].read_header()["EXPOSURE"]
    plate = h[0].read_header()["PLATEID"]
    mjd = h[0].read_header()["MJD"]
    flavor = h[0].read_header()["FLAVOR"]

    flux = (h[0].read() + h[6].read())*h[8].read()*spFlat[0].read()

    ma = h[2].read()
    ## turn off bits coming from sky subtraction:
    ma = ma & (~2**22) & (~2**23) & (~2**27)
    ## turn off bits coming from flat-fielding
    ma = ma & (~2**17) & (~2**2)

    ivar = h[1].read()*(ma==0)*h[8].read()**2*spFlat[0].read()**2
    coeff_wave = h[3]["COEFF"][:]
    coeff_wdisp = h[4]["COEFF"][:]

    if fibers is not None:
        flux = flux[fibers]
        ivar = ivar[fibers]
        coeff_wave = coeff_wave[:,fibers,:]
        coeff_wdisp = coeff_wdisp[:,fibers,:]


    x = np.arange(flux.shape[1])

    spectra = []
    for i in range(flux.shape[0]):
        wave = 10**parse_coeff(x,coeff_wave[0,i,:],h[3][0])
        wdisp = parse_coeff(x,coeff_wdisp[0,i,:],h[4][0])
        spec = spectrum(wave,flux[i],ivar[i],ma[i],wdisp)
        spectra.append(spec)

    return frame(plate, mjd, c, expid, flavor, spectra)

def export_frame(frame,header,fibermap,colnames,prefix):
    
    dirout = prefix+"/exposures/{}/{}".format(header["NIGHT"],frame.expid)
    header["EXPID"] = frame.expid
    header["FLAVOR"] = frame.flavor
    try:
        makedirs(dirout)
    except:
        pass

    spectro = frame.camera[1]

    assert spectro == '1' or spectro == '2'

    if spectro == '1':
        fibermap = [f[:500] for f in fibermap]
    elif spectro == '2':
        fibermap = [f[500:] for f in fibermap]
    
    str_expid = str(frame.expid)
    while len(str_expid)<8:
        str_expid = '0'+str_expid
    fout_name = dirout+"/frame-{}-{}.fits".format(frame.camera,str_expid)
    fout = fitsio.FITS(fout_name,"rw",clobber=True)

    flux = [s.flux for s in frame.spectra]
    ivar = [s.ivar for s in frame.spectra]
    mask = [s.mask for s in frame.spectra]
    reso = [s.wdisp for s in frame.spectra]

    head0 = {}
    head0["CAMERA"] = frame.camera
    head0["FLAVOR"] = frame.flavor

    fout.write(np.array(flux),extname="FLUX",header=head0)
    fout.write(np.array(ivar),extname="IVAR")
    fout.write(np.array(mask).astype(int),extname="MASK")
    fout.write(frame.spectra[0].wave,extname="WAVELENGTH")
    fout.write(np.array(reso),extname="RESOLUTION")
    fout.write(fibermap,header=header,names=colnames,extname="FIBERMAP")

    fout.close()
