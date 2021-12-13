# Troubleshooting 
root -l -b -q 'makeDistributionPlots.C("RelValElectronGunPt2To100_NoPU_vanilla", \
                                       "RelVal Electron Gun Pt 2 to 100, No PU",\
                                       "/Users/stephaniekwan/Dropbox/Princeton_G4/Phase2RCT/analyzer/analyzer.root")'



# root -l -b -q 'runDistributionPlots.C("Nov_10_2021_no_calibration_EGamma_presentation", \
#                                       "RelVal Electron Gun Pt 2 to 100, No PU",\
#                                       "/Users/stephaniekwan/Dropbox/Princeton_G4/Phase2RCT/analyzer/analyzer_11_10_21_no_calibration_EGamma_presentation.root")'

# Regular plots
# root -l -b -q 'runDistributionPlots.C("RelValElectronGunPt2To100_NoPU", \
#                                       "RelVal Electron Gun Pt 2 to 100, No PU",\
#                                       "/Users/stephaniekwan/Dropbox/Princeton_G4/Phase2RCT/analyzer/analyzer.root")'


root -l -b -q 'makeDistributionPlots.C("RelValElectronGunPt2To100_NoPU_old_Cecile_emulator_vanilla", \
                                       "RelVal Electron Gun Pt 2 to 100, No PU",\
                                       "/Users/stephaniekwan/Dropbox/Princeton_G4/Phase2RCT/analyzer/analyzer-l1egammaCrystals.root")'
