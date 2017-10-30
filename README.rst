=========
boss2desi
=========
Code to convert BOSS data to DESI format
----------------------------------------

Example run::

    python bin/convert-frame --frame frame.fits --cframe spCFrame-b1-00121278.fits --wa --out frame-b1-00121278.fits.gz

This will create the ``frame-b1-00121278.fits.gz`` file following the DESI
brick conventions:

* HDU0: flux
* HDU1: ivar
* HDU2: mask
* HDU3: wavelength
* HDU4: resolution
* HDU5: fibermap

The wavelength is determined by HDU3, fiber 250 (or 750) in the ``--wave``
file. Everything else is interpolated to this grid.

The flux is the sum of the flux and sky as read from the frame file
(to allow the user to run a sky subtraction).

The resolution is gaussian with a dispersion given by wdisp. 

The spCFrame is needed to read the wavelengths (and the fibermap,
although this info is also in the spFrame file).
