from scipy.spatial import ConvexHull
import numpy as np
import astropy.units as u

class BBox:
    def  __init__(self):
        self.corners = None
        self.position_angle = None
        self.area = None
        self.length = None
        self.width = None
        self.xcen = None
        self.ycen = None
        
    def __repr__(self):
        return "Bounding box L:{0}, W:{1} oriented at {2}".format(
            self.length, self.width, self.position_angle)

def rotator(x,y,phi):
    xrot = x * np.cos(phi) - y * np.sin(phi)
    yrot = x * np.sin(phi) + y * np.cos(phi)
    return xrot, yrot

def rotbbox(x, y):
    hull = ConvexHull(np.c_[x,y])
    angle = np.zeros(len(hull.vertices))
    area = np.zeros(len(hull.vertices))
    x0 = np.mean(hull.points[:, 0])
    y0 = np.mean(hull.points[:, 1])
    for idx, vert in enumerate(hull.vertices):
        p1 = hull.points[vert]
        p2 = hull.points[hull.vertices[idx-1]]
        angle[idx] = np.arctan((p2[1]-p1[1])/(p2[0]-p1[0]))
    x = hull.points[:, 0] - x0
    y = hull.points[:, 1] - y0
    for idx, phi in enumerate(angle):
        xrot, yrot = rotator(x, y, -phi)
        DeltaX = np.nanmax(xrot) - np.nanmin(xrot)
        DeltaY = np.nanmax(yrot) - np.nanmin(yrot)
        area[idx] = DeltaX * DeltaY
    phi = angle[np.argmin(area)]
    xrot, yrot = rotator(x, y, -phi)
        
    minxrot, maxxrot = np.nanmin(xrot), np.nanmax(xrot)
    minyrot, maxyrot = np.nanmin(yrot), np.nanmax(yrot)
    length = (maxxrot - minxrot)
    width = (maxyrot - minyrot)
    area = length * width
    
    cornrot = np.array([[minxrot, minyrot],
                        [minxrot, maxyrot],
                        [maxxrot, maxyrot],
                        [maxxrot, minyrot]])
    xcorn, ycorn = rotator(cornrot[:, 0],
                           cornrot[:, 1], phi)
    xcorn += x0
    ycorn += y0
    corners = np.c_[xcorn, ycorn]
    box = BBox()
    box.corners = corners
    box.area = area
    box.length = length
    box.width = width
    box.position_angle = phi
    box.xcen = x0
    box.ycen = y0
    
    return box

