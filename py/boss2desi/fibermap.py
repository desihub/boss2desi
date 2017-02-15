import scipy as sp
import fitsio

class boss_fibermap:
    def __init__(self,fibermap,camera,thid_map=None,plate=None):

        ## get boss_fibermap names
        fm_names = fibermap.get_colnames()
        nspec = len(fibermap[fm_names[0]][:])
        ## create a fiber map and add all columns from boss' fibermap
        fm = []
        for i in fm_names:
            fm.append(fibermap[i][:])

        ## add fiber number
        fm_names.append("FIBER")
        fib_num = sp.arange(nspec,dtype=int)
        spectro = int(camera[1])
        if spectro == 2:
            fib_num += 500
        fm.append(fib_num)

        ## add channel
        fm_names.append("CHANNEL")
        fm.append(sp.array([camera]*nspec))

        ## add filter names
        fm_names.append("FILTER")
        fm.append(sp.array([["SDSS_U","SDSS_G","SDSS_R","SDSS_I","SDSS_Z"]]*nspec))

        fm_names.append("SPECTROID")
        fm.append(sp.array([spectro]*nspec))

        ## copy targetid from THINGID
        if thid_map != None:
            fm_names.append("TARGETID")
            targetid = []
            for f in fib_num:
                targetid.append(thid_map[f])
            fm.append(sp.array(targetid))

        if plate is not None:
            fm_names.append("BRICKNAME")
            fm.append(sp.array([plate]*nspec))

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

        self.camera = camera
        self.fm_names = fm_names
        self.fm = fm

    def __radd__(self,other):
        for i,l in enumerate(self.fm):
            self.fm[i] = sp.concatenate((l,other.fm[i]))
    
        return self
    def export(self,fout):

        brick = fitsio.FITS(fout,"rw",clobber=True)
        brick.write(self.fm,names=self.fm_names,extname="FIBERMAP")
        brick.close()
