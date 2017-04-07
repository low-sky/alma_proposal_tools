import astropy.units as u
from astropy.io import fits
from mosaictools import rotbbox
from galaxies import Galaxy

def contourbound(imagename, level=0.5 * u.MJy/u.sr):
    header = fits.getheader(imagename)
    data = fits.getdata(imagename)
