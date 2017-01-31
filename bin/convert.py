from boss2desi import frame2brick
import fitsio
import argparse


if __name__=="__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("--frame",type=str,default=None,required=True,help="frame file")
    parser.add_argument("--cframe",type=str,default=None,required=True,help="cframe file")
    parser.add_argument("--wave",type=str,default=None,required=True,help="frame file defining the wavelength")
    parser.add_argument("--out",type=str,default=None,required=True,help="output file")
    parser.add_argument("--camera",type=str,default=None,required=False,help="the camera (in case it's not found on hdu[0]")
    parser.add_argument("--flat",action="store_true",required=None,help="set if this is a flat file (will read the wavelength from the file")

    args=parser.parse_args()

    frame=fitsio.FITS(args.frame)
    flux = frame[0].read()
    ivar = frame[1].read()
    mask = frame[2].read()
    sky  = None
    try:
        sky  = frame[6].read()
    except:
        pass
    try:
        camera = frame[0].read_header()["CAMERAS"]
    except:
        camera = args.camera

    ## read wdisp fibermap and wavelenghts and wdisp from cframe
    cframe=fitsio.FITS(args.cframe)
    if args.flat:
        wave = 10**frame[3].read()
        wave_wdisp = 10**cframe[3].read()[:,:flux.shape[1]]
    else:
        wave = 10**cframe[3].read()[:,:flux.shape[1]]
        wave_wdisp = wave

    wdisp= cframe[4].read()[:,:flux.shape[1]]
    fibermap = cframe[5]

    ## interpolate into fiber 250 from flat-field wavelengths
    wave_new = fitsio.FITS(args.wave)
    wave_new = 10**wave_new[3].read()[249,:]

    b=frame2brick.brick(flux,ivar,mask,wave,wave_new,wdisp,camera,fibermap,sky=sky,wave_wdisp=wave_wdisp)

    b.export(args.out)


    
