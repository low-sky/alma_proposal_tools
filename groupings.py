from astropy.table import Table, Column
import sklearn.cluster as cluster
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
import scipy.sparse.csgraph as csg
import matplotlib.pyplot as plt

torig = Table.read('califa_sample.csv')
maxmember=5
lowvel = torig['Vlsr']<50000
highvel = torig['Vlsr']>=50000
ubergroup = np.zeros(len(torig))

inttime_18mjy = {1:(28.54)/60,2:39.13/60,3:49.72/60,4:1.01,5:1.18}
inttime_10mjy = {4:1.17}
for offset,thisset in enumerate([lowvel,highvel]):
    t = torig[thisset]
    coords = SkyCoord(ra=t['RA']*15,dec=t['DE'],unit=(u.deg,u.deg))
    d = np.zeros((len(coords),len(coords)))
    dv = t['Vlsr'].data
    dv+= (1,)
    dv = np.abs(dv-dv.T)
    for idx,row in enumerate(coords):
        d[:,idx] = (coords.separation(coords[idx])).value
    d_original = np.copy(d)
    d[d>15.0] = 0.0
    d[dv>3500] = 0.0
    graph=csg.csgraph_from_dense(d)
    ncomp,asgn = csg.connected_components(graph,directed=False)
    group = np.zeros(len(t))


    for thiscomp in np.arange(ncomp):
        subset = (asgn == thiscomp)
        group[subset] = thiscomp*100+offset*1000
        subarray_original = d[np.ix_(subset,subset)]
        subarray = np.copy(subarray_original)
        tsub = t[subset]
        subgroupcounter = 0
        subgroup = np.zeros(subset.sum())
        subgroupidx_orig = np.arange(len(subgroup))
        subgroupidx = np.copy(subgroupidx_orig)

        while np.any(subgroup==0):
            subgroupcounter += 1
            subgraph = csg.csgraph_from_dense(subarray)
            mst = csg.minimum_spanning_tree(subgraph)
            g1,g2 = np.unravel_index(np.argmax(mst.toarray()),subarray.shape)
            treefromhere,predecessors = csg.breadth_first_order(mst,g1,directed=False)
            closestmembers = treefromhere[0:maxmember]
            try:
                dmat = subgraph[np.ix_(closestmembers,closestmembers)]
            except ValueError:
                dmat = subgraph[0]
            reduce = 0
            while np.any(dmat.toarray()>15):
                reduce+=1
                closestmembers = treefromhere[0:(maxmember-reduce)]
                dmat = subgraph[np.ix_(closestmembers,closestmembers)]
            subgroup[subgroupidx[treefromhere[0:(maxmember-reduce)]]] = subgroupcounter
            subgroupidx  = subgroupidx_orig[subgroup == 0] 
            subarray = subarray_original[np.ix_(subgroup == 0,subgroup == 0)]
            print(subgroup)
        group[subset]+=subgroup
    ubergroup[thisset] = group
for relabel,value in enumerate(np.unique(ubergroup)):
    ubergroup[ubergroup==value]=relabel+1
ngal,edges,_ = plt.hist(ubergroup,bins=np.arange(ubergroup.max()+1)+0.5)
totaltime = 0.0
for val in ngal:
    totaltime+=inttime_18mjy[np.int(val)]

nom = torig['Name']
coords = SkyCoord(ra=torig['RA']*15,dec=torig['DE'],unit=(u.deg,u.deg))
cstr = coords.to_string(style='hmsdms',sep=':')
cstr = [c.replace(' ',', ') for c in cstr]
for n, c in zip(nom,cstr):
    print("{0}, {1}, 0.0, 0.0, 5400.0, hel, RELATIVISTIC, 0.0, 0.0, 0.0, 0.0, 0.0".format(n,c))

#links = csg.minimum_spanning_tree(d)

