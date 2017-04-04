import numpy as np

pbwidth = 0.51093*51.013#1.22*(3e8/(115.2712e9*12))

x = np.zeros(7)
y = np.zeros(7)
phi = np.arange(6)*np.pi/3
x[1:] = pbwidth*np.cos(phi)
y[1:] = pbwidth*np.sin(phi)
hdr = 'RA , Dec, Coordinate Type, Coordinate Units\n--   This signals end of the header\n'
ptfile = open("alma_ptfile.txt","w")
ptfile.write(hdr)
for thisx,thisy in zip(x,y):
    ptfile.write("{0},{1},Offset,ARCSECS\n".format(thisx,thisy))
ptfile.close()
