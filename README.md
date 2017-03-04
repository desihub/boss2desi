## Code to convert BOSS data to desi format

example run:

python bin/convert-frame --frame frame.fits --cframe spCFrame-b1-00121278.fits --wa --out frame-b1-00121278.fits.gz

This will create the frame-b1-00121278.fits.gz file following the desi brick conventions:

hdu0: flux
hdu1: ivar
hdu2: mask
hdu3: wavelength
hdu4: resolution
hdu5: fibermap

The wavelength is determined by hdu3, fiber 250 (or 750) in the --wave file. Everything else is interpolated to this grid.

The flux is the sum of the flux and sky as read from the frame file (to allow the user to run a sky subtraction)

The resolution is gaussian with a dispersion given by wdisp. 

The spCFrame is needed to read the wavelengths (and the fibermap, although this info is also in the spFrame file)


