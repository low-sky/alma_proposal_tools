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

def get_band_id(spectralscan):
    start = float(spectralscan.find('.//prj:startFrequency', root.nsmap).text)
    if start < 116:
        band_id = 3
    elif start < 160:
        band_id = 4
    elif start < 211:
        band_id = 5
    elif start < 275:
        band_id = 6
    elif start < 360:
        band_id = 7

    return band_id

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

    #band_id = spectralscan.find('.//prj:requiredReceiverBands', root.nsmap).text[-1]
    band_id = get_band_id(spectralscan)


    with open('singletuning.xml','r') as fh:
        blank_single_tuning_tree = lxml.etree.parse(fh)
        blank_single_tuning = blank_single_tuning_tree.getroot().find('.//prj:ScienceGoal', root.nsmap)

    tuning_ct = 1

    if_sep = if_seps[band_id]

    lastfrq = start + if_sep + 2*increment
    while lastfrq < end:
        tuning = copy.copy(blank_single_tuning)
        bands = tuning.findall('.//prj:ScienceSpectralWindow', root.nsmap)
        frq1 = start + 1*increment/2.
        frq2 = start + 3*increment/2.
        frq3 = start + 1*increment/2. + if_sep
        frq4 = start + 3*increment/2. + if_sep
        start = start + 2*increment
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

        for ii,(bb,frq) in enumerate(zip(bands, (frq1,frq2,frq3,frq4))):
            bb.find('prj:centerFrequency', root.nsmap).text = str(frq)
            bb.find('prj:transitionName', root.nsmap).text = "SpectralScan {0}".format(ii)
            bb.find('prj:groupIndex', root.nsmap).text = str(ii)

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
        band_ids[band] += n_in_band

for sg in root.findall('.//prj:ScienceGoal', root.nsmap):
    LAS = sg.find('.//prj:desiredLargestScale', root.nsmap)
    LAS.text = str(180.0)

    centerFreq = sg.find('.//prj:centerFrequency', root.nsmap).text

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
#for elt in root.findall('.//prj:desiredSensitivity', root.nsmap):
#    elt.text = '0.15'

#for elt in root.findall('.//prj:PerformanceParameters', root.nsmap):
#    elt.attrib['desiredSensitivityFrequencyMeasure'] = 'User'
#    frq = float(elt.find('prj:representativeFrequency', root.nsmap).text)*u.GHz
#    dv = 1.0 * u.km/u.s
#    dnu = (dv / constants.c * frq).to(u.GHz)
#    frqelt = elt.find('prj:desiredSensitivityReferenceFrequencyWidth', root.nsmap)
#    frqelt.text = str(dnu.value)

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

