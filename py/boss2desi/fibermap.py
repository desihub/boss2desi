import numpy as np
import fitsio
from redrock.external.boss import platemjdfiber2targetid
from astropy.time import Time

def desi_fibermap(spplate):
    '''
    returns a desi fibermap given an spPlate fibermap
    '''

    h = fitsio.FITS(spplate)
    boss_head = h[0].read_header()
    boss_fibermap = h[5]


    nspec = len(boss_fibermap["FIBERID"][:])

    ## create a fiber map and add all columns from boss' fibermap

    header={}
    header["TELRA"] = boss_head["RA"]
    header["TELDEC"]= boss_head["DEC"]
    header["FLAVOR"] = 'dark'
    mjd = Time(boss_head["MJD"],format="mjd")
    year = mjd.datetime.year
    month = mjd.datetime.month
    day = mjd.datetime.day
    YYYYMMDD = day + month*100+year*10000
    header["NIGHT"] = "{}".format(YYYYMMDD)
    header["TILE"] = boss_head["PLATEID"]

    fm_names = []
    fm = []
    

    ## gal extinction
    sfd_ebv=h[5]["SFD_EBV"][:]
    fm_names.append("SFD_EBV")
    fm.append(sfd_ebv)

    ## add objtype
    fm_names.append("OBJTYPE")
    obj = [o.strip().replace("SPECTROPHOTO_STD","STD") for o in boss_fibermap["OBJTYPE"][:]]
    fm.append(np.array(obj))

    ## TARGETCAT
    fm_names.append("TARGETCAT")
    fm.append(np.array(["BOSS"]*nspec))

    ## BRICKNAME
    fm_names.append("BRICKNAME")
    fm.append(np.array([str(boss_head["PLATEID"])]*nspec))

    ## TARGETID
    fm_names.append("TARGETID")

    plate = np.array([boss_head["PLATEID"]]*nspec)
    mjd = np.array([boss_head["MJD"]]*nspec)
    fiber = boss_fibermap["FIBERID"][:]
    fm.append(platemjdfiber2targetid(plate,mjd,fiber))

    ## TARGETIBITS
    fm_names.append("DESI_TARGET")
    fm.append(np.array([0]*nspec))

    fm_names.append("BGS_TARGET")
    fm.append(np.array([0]*nspec))

    fm_names.append("MWS_TARGET")
    fm.append(np.array([0]*nspec))

    ## MAGNITUDES
    mag = np.zeros((nspec,5))
    for i,obj in enumerate(boss_fibermap["MAG"][:]):
        ## avoid NA and SKY
        ## sky fibers don't have magnitudes
        if boss_fibermap['OBJTYPE'][:][i].strip() != 'SKY'\
                    and boss_fibermap['OBJTYPE'][:][i].strip() != 'NA':
            ## fix only positive fluxes...
            w = boss_fibermap['CALIBFLUX'][:][i] > 0
            mag[i,w] = 22.5-2.5*np.log10(boss_fibermap['CALIBFLUX'][:][i][w])

    fm_names.append("MAG")
    fm.append(mag)

    ## add filter names
    fm_names.append("FILTER")
    fm.append(np.array([["SDSS_U","SDSS_G","SDSS_R","SDSS_I","SDSS_Z"]]*nspec))

    ## spectroid is (fiberid-1)%500
    fm_names.append("SPECTROID")
    fm.append((boss_fibermap["FIBERID"][:]-1)%500)

    ## positioner. No idea what this is, set to -1
    fm_names.append("POSITIONER")
    fm.append(np.array([-1]*nspec))

    ## fiber = fiberid-1
    fm_names.append("FIBER")
    fm.append(boss_fibermap["FIBERID"][:]-1)

    ## lambda ref
    fm_names.append("LAMBDAREF")
    fm.append(boss_fibermap["LAMBDA_EFF"][:])

    ## RA,DEC
    fm_names.append("RA")
    fm.append(boss_fibermap["RA"][:])
    fm_names.append("DEC")
    fm.append(boss_fibermap["DEC"][:])

    ## RA, DEC deduced from XY, set to RA,DEC
    fm_names.append("RA_OBS")
    fm.append(boss_fibermap["RA"][:])
    fm_names.append("DEC_OBS")
    fm.append(boss_fibermap["DEC"][:])

    ## X_TARGET, Y_TARGET
    fm_names.append("X_TARGET")
    fm.append(boss_fibermap["XFOCAL"][:])
    fm_names.append("Y_TARGET")
    fm.append(boss_fibermap["YFOCAL"][:])

    fm_names.append("X_FVCOBS")
    fm.append(boss_fibermap["XFOCAL"][:])
    fm_names.append("Y_FVCOBS")
    fm.append(boss_fibermap["YFOCAL"][:])

    fm_names.append("X_FVCERR")
    fm.append(boss_fibermap["XFOCAL"][:]*0.01)
    fm_names.append("Y_FVCERR")
    fm.append(boss_fibermap["YFOCAL"][:]*0.01)


    return header,fm_names,fm
