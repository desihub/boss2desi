pro run_extract,yanny_name

;;; read yanny file
print,yanny_name
allseq = yanny_readone(yanny_name,'SPEXP',hdr=hdr,/anon)
mjd=yanny_par(hdr,'MJD')
plate=yanny_par(hdr,'plateid')
w=where(allseq.flavor eq 'arc')
arcnames=allseq[w].name
w=where(allseq.flavor eq 'flat')
flatnames=allseq[w].name
w=where(allseq.flavor eq 'science')
sciencenames=allseq[w].name

;;; find the spArc
spArc = strarr(n_elements(arcnames))
for i=0,n_elements(arcnames)-1 do spArc[i]='spArc'+strsplit(arcnames[i],'sdR',/extract)+'s.gz'

;; find the spFlat names
spFlat = strarr(n_elements(flatnames))
for i=0,n_elements(flatnames)-1 do spFlat[i]='spFlat'+strsplit(flatnames[i],'sdR',/extract)+'s.gz'

for i=0,n_elements(arcnames)-1 do begin
    outname = 'arc'+strsplit(arcnames[i],'sdR',/extract)+'s'
    extract,arcnames[i]+'.gz',spArc[i],spFlat[i],'',outname,mjd=string(mjd),plate=string(plate)
endfor

for i=0,n_elements(flatnames)-1 do begin
    outname = 'flat'+strsplit(flatnames[i],'sdR',/extract)+'s'
    extract,flatnames[i]+'.gz',spArc[i],spFlat[i],'',outname,mjd=string(mjd),plate=string(plate)
endfor

for i=0,n_elements(sciencenames)-1 do begin
    outname = 'frame'+strsplit(sciencenames[i],'sdR',/extract)+'s'
    extract,sciencenames[i]+'.gz',spArc[i mod 4],spFlat[i mod 4],'',outname,mjd=string(mjd),plate=string(plate)
endfor


end
