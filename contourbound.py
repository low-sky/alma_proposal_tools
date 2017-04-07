import astropy.units as u
from astropy.io import fits
from mosaictools import rotbbox
from galaxies import Galaxy
import glob
import matplotlib._cntr as cntr
import matplotlib.path as mplPath
import numpy as np
import matplotlib.pyplot as plt
import astropy.wcs as wcs
from astropy.wcs.utils import proj_plane_pixel_scales as ppps
from astropy.coordinates import SkyCoord
from astropy.table import Table, Column

def imglist(datadir='/mnt/work/erosolow/phangs/proposal/lp_proposal_unwise_grab/'):
    return(glob.glob(datadir + '*band3*fits'))

def contourbound(imagename, level=0.5 * u.MJy / u.sr, doplot=False):
    header = fits.getheader(imagename)
    data = fits.getdata(imagename)
    level = 0.5
    X,Y = np.meshgrid(np.arange(data.shape[0]),
                      np.arange(data.shape[1]))
    c = cntr.Cntr(X, Y, data.T)
    nlist = c.trace(level, level, 0)
    segs = nlist[:len(nlist)//2]
    shortfile = (imagename.split('/'))[-1]
    galname = (shortfile.split('_'))[0]
    g = Galaxy(galname.upper())
    w = wcs.WCS(header)
    xcen, ycen = w.wcs_world2pix(g.center_position.ra.value,
                                 g.center_position.dec.value, 0)
    hits = np.zeros(len(segs), dtype=np.bool)

    for idx, seg in enumerate(segs):
        path = mplPath.Path(seg)
        hits[idx] = path.contains_point((ycen, xcen))
    try:
        seg = segs[np.where(hits)[0][0]]
    except IndexError:
        return(g, None, w)
    box = rotbbox(seg[:, 1], seg[:,0])
    if doplot:
        plt.imshow(data,vmin=0,vmax=5, cmap='Greys', origin='lower')
        plt.plot(seg[:,1], seg[:,0], color='red')
        plt.plot(box.corners[:, 0], box.corners[:, 1])
        plt.show()
    return(g, box, w)
    # import pdb; pdb.set_trace()


#Offset Longitude (arcsec),Offset Latitude (arcsec),p (arcsec),q (arcsec),PositionAngle (deg)
def genboxes():
    ll = imglist()
    t = Table(names=('Galaxy','RA','Dec','LongOff',
                     'LatOff','p','q','PA'),
              dtype=('S8','f8','f8','f8',
                     'f8','f8','f8','f8'))
    for l in ll:
        g, box, www = contourbound(l, doplot=True)
        if box is None:
            print "No contour? ", g.name
        else:    
            t.add_row()
            ra0, dec0 = www.wcs_pix2world(box.xcen, box.ycen, 0)
            mosaicCenter = SkyCoord(ra0, dec0, unit=u.deg, frame='fk5')
            offset = g.center_position.spherical_offsets_to(mosaicCenter)
            dx = ppps(www)
            t[-1]['Galaxy'] = g.name
            t[-1]['RA'] = g.center_position.ra.value
            t[-1]['Dec'] = g.center_position.dec.value
            t[-1]['LongOff'] = offset[0].to(u.arcsec).value
            t[-1]['LatOff'] = offset[1].to(u.arcsec).value
            t[-1]['p'] = box.length * dx[0] * 3600
            t[-1]['q'] = box.width * dx[1] * 3600
            t[-1]['PA'] = box.position_angle * 180 / np.pi
        #        import pdb; pdb.set_trace()
    t.write('boxes.csv')
