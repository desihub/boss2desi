pro extract, inname,arcname,flatname,plugfile,outname,mjd=mjd,plate=plate,final_highrej=final_highrej,final_lowrej=final_lowrej

; need wset xpeak lambda xsol fflat fibermask
; [transpose(lambda), xpeak], arcinfofile


run2d='v5_9_0'
print,"INFO reading arc ",arcname
arcname=getenv('BOSS_SPECTRO_DATA')+'/../redux/'+run2d+'/'+plate+'/'+arcname
print,arcname
toto=mrdfits(arcname,1)
lambda=transpose(toto[0,*])
xpeak=toto[1:*,*]
wset=mrdfits(arcname,2)
traceset2xy, wset, ypixw, loglam

indir=getenv('BOSS_SPECTRO_DATA')+'/'+mjd

print,"INFO reading flat ",flatname
flatname=getenv('BOSS_SPECTRO_DATA')+'/../redux/'+run2d+'/'+plate+'/'+flatname
print,flatname
fflat=mrdfits(flatname,0)
xset=mrdfits(flatname,1)
traceset2xy, xset, ypix, xtrace

fibermask=mrdfits(flatname,2)
widthset=mrdfits(flatname,3)

print,"INFO preprocessing ",inname
sdssproc, inname, image, invvar, indir=indir, hdr=objhdr, rdnoiseimg=rdnoise, $
       spectrographid=spectrographid, color=color, $
       ecalibfile=ecalibfile, minflat=0.8, maxflat=1.2, $
       nsatrow=nsatrow, fbadpix=fbadpix, /applycrosstalk, ccdmask=ccdmask


print,"INFO reading plugmap ",plugfile
plugmap = readplugmap(plugfile, spectrographid, $
    plugdir=plugdir, /calibobj, mjd=sxpar(objhdr,'MJD'), indir=indir, $
    exptime=sxpar(objhdr,'EXPTIME'), hdr=hdrplug)

; copied from extract_object (cannot call because of hardcoded sky line fit)
configuration=obj_new('configuration', sxpar(objhdr, 'MJD'))
objname = strtrim(sxpar(objhdr,'OBJFILE'),2) 
flavor  = strtrim(sxpar(objhdr,'FLAVOR'),2) 
camera  = strtrim(sxpar(objhdr,'CAMERAS'),2) 
fextract = extract_boxcar(image*(invvar GT 0), xtrace)
scrunch = djs_median(fextract, 1) ; Find median counts/row in each fiber
whopping = find_whopping(scrunch, 10000.0, whopct)
scrunch_sort = sort(scrunch)
i5 = n_elements(scrunch)/20
i95 = i5 * 19
splog, 'Whopping fibers: ', whopping
splog, 'Median counts in all fibers = ', djs_median(scrunch)
splog, 'Number of bright fibers = ', whopct

badcheck = extract_boxcar((invvar LE 0), xtrace, radius=2.5)
badplace = where(badcheck GT 0)
nx = (size(fextract,/dim))[0] 
ny = (size(fextract,/dim))[1] 
pixelmask = lonarr(nx,ny)

badcolumns = where(total(badcheck GT 0,1) GT 0.45 * nx) ; change from 0.1 ???

if (badplace[0] NE -1) then pixelmask[badplace] = $
  pixelmask[badplace] OR pixelmask_bits('NEARBADPIXEL')

if (badcolumns[0] NE -1) then fibermask[badcolumns] = $
  fibermask[badcolumns] OR pixelmask_bits('MANYBADCOLUMNS')

if (whopping[0] NE -1) then begin
                                ; Set the mask bit for whopping fibers themselves
    fibermask[whopping] = fibermask[whopping] OR pixelmask_bits('WHOPPER')
    
                                ; Set the mask bit for fibers near whopping fibers, excluding the
                                ; whopping fibers themselves.  Note that a fiber could still have both
                                ; WHOPPER and NEARWHOPPER set if it is both bright and near another
                                ; bright fiber.
    wp = [whopping - 2 , whopping -1, whopping+1 , whopping+2]
    wp = wp[ where(wp GE 0 AND wp LT ny) ]
    fibermask[wp] = fibermask[wp] OR pixelmask_bits('NEARWHOPPER')
endif

                                ;-----
                                ; Inherit any mask bits from the ccdmask, by setting the pixmask bits
                                ; for anything that would be hit in a boxcar extraction

if (keyword_set(ccdmask)) then begin
    for ibit=0, 31 do begin
        thischeck = extract_boxcar((ccdmask AND 2L^ibit) NE 0, xtrace, $
                                   radius=2.5)
        pixelmask = pixelmask OR (2L^ibit * (thischeck GT 0))
    endfor
endif

                                ;-----------------------------------------------------------------------
                                ;  This is a kludge to fix first and last column ???
                                ;-----------------------------------------------------------------------
if (configuration->extract_object_fixcolumns()) then begin
    image[0,*] = image[0,*]*0.7
    image[nx-1,*] = image[nx-1,*]*0.7
endif

                                ;
                                ;  First we should attempt to shift trace to object flexure
                                ;

xnow = match_trace(image, invvar, xtrace)
bestlag = median(xnow-xtrace)

splog, 'Shifting traces by match_trace ', bestlag

if (abs(bestlag) GT 1.0) then begin
    splog, 'WARNING: pixel shift is large!'
endif

highrej = 10                    ; just for first extraction steps
lowrej = 10                     ; just for first extraction steps
                                ; We need to check npoly with new scattered light backgrounds
npoly = 16                   ; maybe more structure, lots of structure
nrow = (size(image))[2]
yrow = lindgen(nrow) 
nfirst = n_elements(yrow)

splog, 'Extracting frame '+objname+' with 4 step process'

traceset2xy, widthset, xx, sigma2
ntrace = (size(sigma2,/dimens))[1]
wfixed = [1,1]   ; Fit gaussian height + width (fixed center position)
nterms = n_elements(wfixed)

                                ;-----------------------------------------------------------------------
                                ;  Now, subtract halo image and do final extraction with all rows
                                ;-----------------------------------------------------------------------
                                ; (6) Final extraction
splog, 'Step 6: Final Object extraction'

if keyword_set(final_highrej) then highrej = final_highrej else highrej = 4
if keyword_set(final_lowrej) then lowrej  = final_lowrej else lowrej=4

wfixed = [1]               ; Fit to height only (fixed width + center)
nterms = n_elements(wfixed)
reject = [0.2,0.2,0.99]
npoly = 0

; ASB: switching to bundle-wise extraction:

print,"INFO size xnow", size(xnow,/dimens)
print,"INFO size sigma2", size(sigma2,/dimens)

extract_bundle_image, image, invvar, rdnoise, xnow, sigma2, flux, fluxivar,$
  proftype=proftype, wfixed=wfixed, ansimage=ansimage3, $
  highrej=highrej, lowrej=lowrej, npoly=2L, $ ; whopping=whopping, $
  chisq=chisq, ymodel=ymodel, pixelmask=pixelmask, reject=reject, /relative,$
  nperbun=20L, buffsize=8L , chi2pdf=chi2pdf


; write this 
outdir=getenv('BOSS_SPECTRO_REDUX')+'/'+run2d+'/'+plate
spawn,'mkdir -p '+outdir
outname=outdir+'/'+outname

sxaddpar, bighdr, 'BUNIT', 'electrons/row'
mwrfits, flux, outname, bighdr, /create
sxaddpar, hdrfloat, 'BUNIT', 'electrons/row'
sxaddpar, hdrfloat, 'EXTNAME', 'IVAR', ' Inverse variance'
mwrfits, fluxivar, outname, hdrfloat
sxaddpar, hdrfloat, 'EXTNAME', 'MASK', ' Inverse variance'
mwrfits, pixelmask, outname, hdrfloat
sxaddpar, hdrfloat, 'BUNIT', 'log10(Angs)'
sxaddpar, hdrfloat, 'EXTNAME', 'LOGLAM', ' log wavelength'
mwrfits, loglam, outname, hdrfloat
sxaddpar, hdrfloat, 'BUNIT', 'NODIM'
sxaddpar, hdrfloat, 'EXTNAME', 'CHI2PDF', ' chi2/ndf'
mwrfits, chi2pdf, outname, hdrfloat



; called in spreduce.pro
; extract_object, outname, objhdr, image, invvar, plugmap, wset, $
;   xpeak, lambda, xsol, fflat, fibermask, color=color, $
;   proftype=proftype, superflatset=superflatset, $
;   widthset=widthset, dispset=dispset, skylinefile=fullskyfile, $
;   plottitle=plottitle, do_telluric=do_telluric, bbspec=bbspec, $
;   splitsky=splitsky, ccdmask=ccdmask


end
