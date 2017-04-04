from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np

t = Table.read('s.csv')
vmid = np.median(t['Vlsr'])
lowv = t[t['Vlsr']<vmid]
highv = t[t['Vlsr']>=vmid]

lowvtune = 0.5*(lowv['Vlsr'].min()+lowv['Vlsr'].max())
highvtune = 0.5*(highv['Vlsr'].min()+highv['Vlsr'].max())

hdr = ["Name, RA(sex), Dec(sex), PMRA(mas/yr), PMDec(mas/yr), vel(km/s), ref frame, Doppler type, peak cont flux(mJy), peak line flux(mJy), cont pol(%), line pol(%), line width(km/s)\n","--   This signals end of the header\n"]

thisline = "0.0, 0.0, {0:.4g}, lsrk, RADIO, 0.0, 0.0, 0.0, 0.0, 0.0\n".format(lowvtune)
lowfile = open("alma_lowcat.txt","w")
for hdrline in hdr:
    lowfile.write(hdrline)
for thisobj in lowv:
    coords = SkyCoord(ra=thisobj['RA']*15,dec=thisobj['DE'],unit=(u.deg,u.deg))
    coordstring = coords.to_string(style='hmsdms',sep=':')
    coordstring = coordstring.replace(" ",", ")
    name = thisobj['Name']
    name = name.replace("-","_")
    thisstring = "{0}, {1}, ".format(name,coordstring)+thisline
    lowfile.write(thisstring)
lowfile.close()


thisline = "0.0, 0.0, {0:.4g}, lsrk, RADIO, 0.0, 0.0, 0.0, 0.0, 0.0\n".format(highvtune)
highfile = open("alma_highcat.txt","w")
for hdrline in hdr:
    highfile.write(hdrline)
for thisobj in highv:
    coords = SkyCoord(ra=thisobj['RA']*15,dec=thisobj['DE'],unit=(u.deg,u.deg))
    coordstring = coords.to_string(style='hmsdms',sep=':')
    coordstring = coordstring.replace(" ",", ")
    name = thisobj['Name']
    name = name.replace("-","_")
    thisstring = "{0}, {1}, ".format(name,coordstring)+thisline
    highfile.write(thisstring)
highfile.close()
