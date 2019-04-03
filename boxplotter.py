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
import aplpy

def imglist(datadir='/mnt/work/erosolow/phangs/proposal/lp_proposal_unwise_grab/'):
    return(glob.glob(datadir + '*band3*fits'))

def dsslist(datadir='/mnt/work/erosolow/phangs/proposal/lp_proposal_unwise_grab/'):
    return(glob.glob(datadir + '*dss*fits'))


def boxplot():
    t = Table.read('boxes_populated.csv')
    
    wise3 = imglist()
    dss = dsslist()
    for row in t:
        name = row['Galaxy']
        wise3name = ([l for l in wise3 if name in l])[0]
        dssname = ([l for l in dss if name in l])[0]
        xcorner = np.array([0.5,0.5,-0.5, -0.5, 0.5]) * row['p'] / 3600
        ycorner = np.array([0.5,-0.5, -0.5, 0.5, 0.5]) * row['q'] / 3600
        xnew, ynew = rotator(xcorner, ycorner,
                             -row['PA'] * np.pi/180)
        cosdec = np.cos(row['Dec'] * np.pi/180)
        xnew /= cosdec
        xnew += (row['LongOff']/3600./cosdec + row['RA'])
        ynew += (row['LatOff']/3600. + row['Dec'])
        fig = plt.figure(figsize=(8,3.5))
        f = aplpy.FITSFigure(dssname, figure=fig,
                             subplot = [0.1, 0.15,0.4,0.73])
        f.show_grayscale(stretch='arcsinh')
        f.show_polygons([np.c_[xnew, ynew]], color='red', lw=2)
        f.set_tick_labels_font(size='x-small')
        f.set_axis_labels_font(size='small')
        g = aplpy.FITSFigure(wise3name, figure=fig,
                             subplot = [0.55, 0.15,0.4,0.73])
        g.show_grayscale(stretch='arcsinh')
        g.set_tick_labels_font(size='x-small')
        g.set_axis_labels_font(size='small')
        g.show_polygons([np.c_[xnew, ynew]], color='red', lw=2)
        fig.suptitle(name+' Subblock {0}'.format(row['SB']))
        fig.savefig(name+'_'+str(row['SB'])+'.png')


