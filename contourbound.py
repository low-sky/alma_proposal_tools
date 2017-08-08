import astropy.units as u
from astropy.io import fits
from mosaictools import rotbbox, rotator, BBox
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

def contourbound(imagename, level=0.5 * u.MJy / u.sr,
                 doplot=False, PixelMask=False):
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
    try:
        g = Galaxy(galname.upper())
    except:
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
        plt.imshow(np.sqrt(data),vmin=0,vmax=5**0.5,
                   cmap='Greys', origin='lower')
        plt.axis([data.shape[0]*0.2, data.shape[0]*0.8,
                  data.shape[1]*0.2, data.shape[1]*0.8])
        plt.plot(seg[:,1], seg[:,0], color='red')
        plt.plot(np.r_[box.corners[:, 0],box.corners[0,0]],
                 np.r_[box.corners[:, 1],box.corners[0,1]])
        plt.savefig(g.name+'_mosaic.png')
        plt.close()
        plt.clf()
    if PixelMask:
        mask = np.zeros_like(data, dtype='bool')
        thispath = mplPath.Path(box.corners)
        mask[X.ravel(), Y.ravel()] = thispath.contains_points(np.c_[X.ravel(),Y.ravel()])
        mask = mask.astype(np.int)
        hdu = fits.PrimaryHDU(mask, header=w.to_header())
        hdu.writeto(imagename.replace('band3','band3_obsmask'), clobber=True)
    return(g, box, w)
    # import pdb; pdb.set_trace()


#Offset Longitude (arcsec),Offset Latitude (arcsec),p (arcsec),q (arcsec),PositionAngle (deg)

def genboxes():
    PixelMask = False
    ll = imglist()
    namelist = np.loadtxt('proposed_targets_v2.txt',dtype='str')
    if len(namelist.shape) > 1:
        namelist = namelist[:, 0]
    t = Table(names=('Galaxy','SB','RA','Dec','LongOff',
                     'LatOff','p','q','PA'),
              dtype=('S8','i4','f8','f8','f8',
                     'f8','f8','f8','f8'))
    for name in namelist:
        l = [l for l in ll if name in l]
        g, box, www = contourbound(l[0], doplot=True)
        if box is None:
            print "No contour? ", g.name
            box = BBox()
            box.corners = None
            box.position_angle = 0.0 * u.deg
            box.area = 3600.
            box.length = 50.0
            box.width = 50.0
            box.xcen = 0.0
            box.ycen = 0.0

        else:    
            ra0, dec0 = www.wcs_pix2world(box.xcen, box.ycen, 0)
            mosaicCenter = SkyCoord(ra0, dec0, unit=u.deg, frame='fk5')
            offset = g.center_position.spherical_offsets_to(mosaicCenter)
            dx = ppps(www)
            p = box.length * dx[0] * 3600
            q = box.width * dx[0] * 3600
            area = p * q
            ntile = np.ceil(p * q / 1.8e4)
            if ntile == 1:
                t.add_row()
                t[-1]['Galaxy'] = g.name
                t[-1]['SB'] = 1
                t[-1]['RA'] = g.center_position.ra.value
                t[-1]['Dec'] = g.center_position.dec.value
                t[-1]['LongOff'] = offset[0].to(u.arcsec).value
                t[-1]['LatOff'] = offset[1].to(u.arcsec).value
                t[-1]['p'] = p
                t[-1]['q'] = q
                t[-1]['PA'] = box.position_angle * 180 / np.pi
                if PixelMask:
                    writepixelmask(l[0], www, t[-1])
#                if 'IC1954' in g.name:
#                    import pdb; pdb.set_trace()
            if ntile > 1:
                PisBigger = p > q
                if PisBigger:
                    Pnew = p / ntile
                    Qnew = q
                else:
                    Pnew = q / ntile
                    Qnew = p
                    box.position_angle += np.pi / 2
                
                Poffset = np.linspace(-(ntile - 1) / 2,
                                      (ntile - 1) / 2, ntile) * Pnew
                Qoffset = np.zeros_like(Poffset)
                x0 = offset[0].to(u.arcsec).value
                y0 = offset[1].to(u.arcsec).value
                xnew, ynew = rotator(Poffset, Qoffset,
                                     box.position_angle)
                xnew = x0 - xnew
                ynew += y0
                for i in np.arange(ntile):
                    t.add_row()
                    t[-1]['Galaxy'] = g.name
                    t[-1]['SB'] = i + 1
                    t[-1]['RA'] = g.center_position.ra.value
                    t[-1]['Dec'] = g.center_position.dec.value
                    t[-1]['LongOff'] = xnew[int(i)]
                    t[-1]['LatOff'] = ynew[int(i)]
                    t[-1]['p'] = Pnew
                    t[-1]['q'] = Qnew
                    t[-1]['PA'] = box.position_angle * 180 / np.pi
                    if PixelMask:
                        writepixelmask(l[0], www, t[-1])
        t.write('boxes.csv')
