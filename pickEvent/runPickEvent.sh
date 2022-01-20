# source runPickEvent.sh
# Takes an input file, the run:lumi:eventnumber to select, and outputs a file with just that event
# First make an n-tuple to decide which run/lumi/event numbers you want to keep, then use this for 
# further studies.

# From https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookDataSamples#Introduction_to_copy_PickMerge_c
#      "How to copy a particular event"

cmsRun pickEvent_cfg.py inputFiles=file:/afs/cern.ch/work/s/skkwan/public/phase2RCT/RelValElectronGunPt2To100_190EDE9F-770B-174A-8BA6-F7814FC67FD4.root eventsToProcess=1:83:8257 outputFile=file:/afs/cern.ch/work/s/skkwan/public/phase2RCT/RelValElectronGunPt2To100_190EDE9F-770B-174A-8BA6-F7814FC67FD4_evt_1_83_8257.root

