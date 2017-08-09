import os
import xml.etree.ElementTree as ET
import subprocess
from astropy.table import Table
import numpy as np

filein = "loaded.aot"
fileout = "loaded_hack.aot"

t = Table.read('s.csv')

nsdict = {'ns0':"Alma/ObsPrep/ObsProposal",
          'ns2':"Alma/ObsPrep/ObsProject",
          'ns3':'Alma/ValueTypes'}

subprocess.call("unzip -o {0}".format(filein),shell=True)

tree = ET.parse('ObsProposal.xml')
root = tree.getroot()
for thisgal in root.findall('.//ns2:TargetParameters',nsdict):
    name = thisgal.find('.//ns2:sourceName',nsdict).text
    print(name)
    catalog_entry = (t[(t['Name'])==name.replace('_','-')])[0]
    pkcflux = thisgal.find('.//ns2:expectedPeakFluxDensity',nsdict)
    pkcflux.text = "0.1"
    pkflux = thisgal.find('.//ns2:expectedPeakLineFluxDensity',nsdict)
    if np.isfinite(catalog_entry['PeakDetJykms']/50.0):
        pkflux.text = "{0:.0f}".format(catalog_entry['PeakDetJykms']/50.0*1e3)
        print(pkflux.text)
    else:
        pkflux.text = "45.0"
    reffreq = np.float(thisgal.find('.//ns2:referenceFrequency',nsdict).text)
    
    linewidth = thisgal.find('.//ns2:expectedLineWidth',nsdict)
    linewidth.text = "{0}".format(50.0/3e5*reffreq)
    palong = thisgal.find('.//ns2:pALong',nsdict)
    palong.text = "{0:.0f}".format(catalog_entry['PA']-90)
    short = thisgal.find('.//ns2:short',nsdict)
    short.text = '67.0'
    lng = thisgal.find('.//ns2:long',nsdict)
    lng.text = '67.0'

    
#    import pdb; pdb.set_trace()

    
#    if 'ObsProject}long' in n.tag:
#        n.text = '60.0'
#    if 'ObsProject}short' in n.tag:
#        n.text = '60.0'

tree.write('ObsProposal.xml')
subprocess.call("zip {0} ObsProject.xml ObsProposal.xml".format(fileout),shell=True)
