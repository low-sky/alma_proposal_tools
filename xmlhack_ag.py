import os
import copy
#import xml.etree.ElementTree as ET
import lxml.etree
from astropy.table import Table
import numpy as np
import zipfile

from astropy import units as u, constants



#def next_band(sciencegoal, idnum, increment=1.875*2):
#    for representativefreq in sciencegoal.findall('.//prj:representativeFrequency', root.nsmap):
#        representativefreq.text = str(float(representativefreq.text) + increment)
#
#    for centerfreq in sciencegoal.findall('.//prj:centerFrequency', root.nsmap):
#        centerfreq.text = str(float(centerfreq.text) + increment)
#
#    name = sciencegoal.find('prj:name', root.nsmap)
#    name.text = 'B3 tuning {0} in G0.253'.format(idnum)
#
#    return sciencegoal

frq_covg = open('frequency_coverage.txt', 'w')

band_edges = {3: (84, 116),
              4: (125, 163),
              5: (163, 211),
              6: (211, 275),
              7: (275, 373),
             }

def get_band_id(spectralscan):
    start = float(spectralscan.find('.//prj:startFrequency', root.nsmap).text)
    for band_id,(lower,upper) in band_edges.items():
        if (start < upper) and (start >= lower):
            return band_id
    raise ValueError("No valid band found")

if_seps = {3: 12,
           4: 12,
           5: 12,
           6: 14,
           7: 12,
          }

def SpectralScanToIndividualSbs(spectralscan, increment=1.875,
                                n_in_band=0):

    root = spectralscan.getroottree().getroot()

    start = float(spectralscan.find('.//prj:startFrequency', root.nsmap).text) + increment / 2.
    end = float(spectralscan.find('.//prj:endFrequency', root.nsmap).text)
    print("Converting spectral scan that covers {0} to {1} GHz".format(start,end))

    #band_id = spectralscan.find('.//prj:requiredReceiverBands', root.nsmap).text[-1]
    band_id = get_band_id(spectralscan)


    with open('singletuning.xml','r') as fh:
        blank_single_tuning_tree = lxml.etree.parse(fh)
        blank_single_tuning = blank_single_tuning_tree.getroot().find('.//prj:ScienceGoal', root.nsmap)

    tuning_ct = 1

    if_sep = if_seps[band_id]

    if (end-start) <= (if_sep + 2*increment):
        lastfrq = start + 2*increment
    else:
        lastfrq = start + if_sep + 2*increment
    assert lastfrq < end
    while lastfrq < end:
        tuning = copy.copy(blank_single_tuning)
        bands = tuning.findall('.//prj:ScienceSpectralWindow', root.nsmap)
        frq1 = start + 1*increment/2.
        frq2 = start + 3*increment/2.
        frq3 = start + 1*increment/2. + if_sep
        frq4 = start + 3*increment/2. + if_sep
        start = start + 2*increment
        if (end-start) <= (if_sep + 2*increment):
            lastfrq = start + 2*increment
        else:
            lastfrq = start + if_sep + 2*increment

        print("start: {0} end: {1} band: {2} start: {3}"
              "    frqs={4}   frq4-frq1={5}  frq3-frq2={6}"
              .format(start, end, band_id, n_in_band,
                      (frq1,frq2,frq3,frq4),
                      frq4-frq1, frq3-frq2))

        for reffrq in tuning.findall('.//prj:referenceFrequency', root.nsmap):
            reffrq.text = str(frq2)
        for repfrq in tuning.findall('.//prj:representativeFrequency', root.nsmap):
            repfrq.text = str(frq2)
        tuning.find('.//prj:requiredReceiverBands', root.nsmap).text = str("ALMA_RB_0{0}".format(band_id))

        if (frq4 > band_edges[band_id][0]) and (frq4 < band_edges[band_id][1]):
            frqs = (frq1, frq2, frq3, frq4)
            bandwidthGHz = 2
            bw_incr = increment
        else:
            # these are duplicate coverage
            continue
            frqs = (frq1, frq1+increment/2., frq2, frq2+increment/2.)
            bandwidthGHz = 1
            bw_incr = increment/2

        for ii,(bb,frq) in enumerate(zip(bands, frqs)):
            bb.find('prj:centerFrequency', root.nsmap).text = str(frq)
            bb.find('prj:transitionName', root.nsmap).text = "SpectralScan {0}".format(ii)
            bb.find('prj:groupIndex', root.nsmap).text = str(ii)

            repwin = bb.find('prj:representativeWindow', root.nsmap)
            if ii == 0:
                repwin.text = 'true'
            else:
                repwin.text = 'false'

            bandwidth = bb.find('prj:bandWidth', root.nsmap)
            bandwidth.text = str(bandwidthGHz)
            bandwidth.attrib['unit'] = 'GHz'
            specres = bb.find('prj:spectralResolution', root.nsmap)
            specres.text = str(bandwidthGHz/4096.)
            specres.attrib['unit'] = 'GHz'

            frq_covg.write("{0} {1}\n".format(frq-bw_incr/2., frq+bw_incr/2.))

        angres = tuning.find('.//prj:desiredAngularResolution', root.nsmap)
        if frq1 < 340:
            angres.text = '1.0'
        else:
            angres.text = '0.75'

        name = tuning.find('prj:name', root.nsmap)
        name.text = 'B{0} tuning {1} G0.253'.format(band_id, tuning_ct+n_in_band)

        tuning_ct += 1

        spectralscan.addnext(tuning)

    root.remove(spectralscan)



    return band_id, n_in_band+tuning_ct-1



#tree = ET.parse('ObsProposal.xml')
#root = tree.getroot()

with zipfile.ZipFile('cmz_linesurvey_lowspectralres_scans.aot','r') as aot:
    aot.extractall()
    os.remove('ObsProposal.xml')
    aot.extract('ObsProject.xml')
    aot.extract('ObsProposal.xml')
    files = aot.filelist

with open('ObsProposal.xml','r') as fh:
    tree = lxml.etree.parse(fh)
    root = tree.getroot()


band_ids = {3: 0, 4: 0, 5: 0, 6: 0, 7: 0}


for ii,spectralscan in enumerate(root.findall('.//prj:ScienceGoal', root.nsmap)):
    print(ii)
    if spectralscan.find('.//prj:SpectralScan', root.nsmap):
        band = get_band_id(spectralscan)
        band_id, n_in_band = SpectralScanToIndividualSbs(spectralscan,
                                                         n_in_band=band_ids[band])
        band_ids[band] = n_in_band
    elif spectralscan.attrib['mode'] == 'Standard':
        for band in spectralscan.findall('.//prj:ScienceSpectralWindow', root.nsmap):

            frqobj = band.find('prj:centerFrequency', root.nsmap)
            frq = u.Quantity(float(frqobj.text), frqobj.attrib['unit']).to(u.GHz).value
            bwobj = band.find('prj:bandWidth', root.nsmap)
            bw_incr = u.Quantity(float(bwobj.text)*1.875/2.0, bwobj.attrib['unit']).to(u.GHz).value
            print("Adding to frequency coverage: ",frq,bw_incr)

            frq_covg.write("{0} {1}\n".format(frq-bw_incr/2., frq+bw_incr/2.))
    else:
        print("Science goal {0} is neither SpectralScan nor normal".format(spectralscan))

for sg in root.findall('.//prj:ScienceGoal', root.nsmap):
    LAS = sg.find('.//prj:desiredLargestScale', root.nsmap)
    LAS.text = str(180.0)

    repwindow = [x for x in sg.findall('.//prj:ScienceSpectralWindow', root.nsmap)
                 if x.find('prj:representativeWindow', root.nsmap).text == 'true'
                ][0]
    centerFreq = repwindow.find('.//prj:centerFrequency', root.nsmap).text

    for reffrq in sg.findall('.//prj:referenceFrequency', root.nsmap):
        reffrq.text = centerFreq
    for repfrq in sg.findall('.//prj:representativeFrequency', root.nsmap):
        repfrq.text = centerFreq

#b3tuning = root.findall('.//prj:ScienceGoal', root.nsmap)[-1]
#assert 'B3 tuning 1' in b3tuning.find('prj:name', root.nsmap).text
#for ii in range(2,7):
#    prev = b3tuning
#    b3tuning = copy.copy(b3tuning)
#    next_band(b3tuning, ii)
#    prev.addnext(b3tuning)


#for elt in root.findall('.//prj:smoothingFactor', root.nsmap):
#    elt.text = '1'
#
for elt in root.findall('.//prj:desiredSensitivity', root.nsmap):
    elt.text = '0.25'

for elt in root.findall('.//prj:PerformanceParameters', root.nsmap):
    elt.attrib['desiredSensitivityFrequencyMeasure'] = 'User'
    frq = float(elt.find('prj:representativeFrequency', root.nsmap).text)*u.GHz
    dv = 1.0 * u.km/u.s
    dnu = (dv / constants.c * frq).to(u.GHz)
    frqelt = elt.find('prj:desiredSensitivityReferenceFrequencyWidth', root.nsmap)
    frqelt.text = str(dnu.value)

#for elt in root.findall('.//prj:SpectralScan', root.nsmap):
#    bandwidth = elt.find('prj:bandWidth', root.nsmap)
#    bandwidth.text = '2.0'
#    specres = elt.find('prj:spectralResolution', root.nsmap)
#    specres.text = str(2.0/4096.)

with open('ObsProposal.xml','w') as fh:
    tree.write('ObsProposal.xml')

with zipfile.ZipFile('cmz_linesurvey_lowspectralres_individualSBs.aot','w') as aot:
    for fh in files:
        aot.write(fh.filename)

frq_covg.close()
