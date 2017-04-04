import os
import xml.etree.ElementTree as ET
import subprocess

filein = "sevenpoint.aot"
fileout = "sevenpoint_hack.aot"
ns = {'prj': 'Alma/ObsPrep/ObsProject', 'prp': 'Alma/ObsPrep/ObsProposal'}

subprocess.call("unzip -o {0}".format(filein),shell=True)
tree = ET.parse('ObsProposal.xml')
root = tree.getroot()
# for child in root:
#     if 'ScienceGoal' in child.tag:
#         for grandchild in child:
#             if 'TargetParameters' in grandchild.tag:
#                 for ggchild in grandchild:
#                     if 'SinglePoint' in ggchild.tag:
#                         for gggchild in ggchild:
#                             print (gggchild.tag)
#                             import pdb; pdb.set_trace()

                    



#tree.write('ObsProposal.xml')
#subprocess.call("zip {0} ObsProject.xml ObsProposal.xml".format(fileout),shell=True)
