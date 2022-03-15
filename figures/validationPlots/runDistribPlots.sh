# runDistribPlot.sh
# Usage: source runDistribPlot.sh


### The usual reco/gen pT cuts, but with different scale factors
root -l -b -q 'makeDistributionPlots.C("SF_1.0_cPt-gt-15_genPt-gt-15",                                               \
                                       "RelVal Electron Gun Pt 2 to 100",                                     \
                                       "/Users/stephaniekwan/Dropbox/Princeton_G4/Phase2RCT/analyzer/analyzer.root", \
                                       "(gct_cPt > 15) && (genPt > 15)", "SF = 1.0",                                 \
                                       "L1 p_{T} > 15, Gen p_{T} > 15",                                              \
                                       "1.0")' 


root -l -b -q 'makeDistributionPlots.C("SF_0.985_cPt-gt-15_genPt-gt-15",                                             \
                                       "RelVal Electron Gun Pt 2 to 100",                                     \
                                       "/Users/stephaniekwan/Dropbox/Princeton_G4/Phase2RCT/analyzer/analyzer.root", \
                                       "(gct_cPt > 15) && (genPt > 15)", "SF = 0.985",                               \
                                       "L1 p_{T} > 15, Gen p_{T} > 15",                                              \
                                       "0.985")'

root -l -b -q 'makeDistributionPlots.C("SF_0.963_cPt-gt-15_genPt-gt-15",                                             \
                                       "RelVal Electron Gun Pt 2 to 100",                                     \
                                       "/Users/stephaniekwan/Dropbox/Princeton_G4/Phase2RCT/analyzer/analyzer.root", \
                                       "(gct_cPt > 15) && (genPt > 15)", "SF = 0.963",                               \
                                       "L1 p_{T} > 15, Gen p_{T} > 15",                                              \
                                       "0.963")' 

### Stick to SF = 0.985, try different pT cuts
root -l -b -q 'makeDistributionPlots.C("SF_0.985_cPt-gt-20_genPt-gt-30",                                             \
                                       "RelVal Electron Gun Pt 2 to 100",                                     \
                                       "/Users/stephaniekwan/Dropbox/Princeton_G4/Phase2RCT/analyzer/analyzer.root", \
                                       "(gct_cPt > 20) && (genPt > 30)", "SF = 0.985",                               \
                                       "L1 p_{T} > 20, Gen p_{T} > 30",                                              \
                                       "0.985")'

root -l -b -q 'makeDistributionPlots.C("SF_0.985_cPt-gt-20-lt-75_genPt-gt-30-lt70",                           \
                                       "RelVal Electron Gun Pt 2 to 100",                                     \
                                       "/Users/stephaniekwan/Dropbox/Princeton_G4/Phase2RCT/analyzer/analyzer.root", \
                                       "(gct_cPt > 20) && (gct_cPt < 75) && (genPt > 30) && (genPt < 70)", "SF = 0.985",                               \
                                       "20 < L1 p_{T} < 75, 30 < Gen p_{T} < 70",                                              \
                                       "0.985")'                                       

root -l -b -q 'makeDistributionPlots.C("SF_0.985_cPt-gt-75_genPt-gt-70",                                      \
                                       "RelVal Electron Gun Pt 2 to 100",                                     \
                                       "/Users/stephaniekwan/Dropbox/Princeton_G4/Phase2RCT/analyzer/analyzer.root", \
                                       "(gct_cPt > 75) && (genPt > 70)", "SF = 0.985",                               \
                                       "L1 p_{T} > 75, Gen p_{T} > 70",                                              \
                                       "0.985")'                                       


# ### Stick to SF = 0.985, just separate good agreement from bad agreement
# root -l -b -q 'makeDistributionPlots.C("SF_0.985_cPt-gt-15_genPt-gt-15_resolution_worse_than_-0.05",                 \
#                                        "RelVal Electron Gun Pt 2 to 100",                                     \
#                                        "/Users/stephaniekwan/Dropbox/Princeton_G4/Phase2RCT/analyzer/analyzer.root", \
#                                        "(gct_cPt > 15) && (genPt > 15) && ((gct_cPt/0.985 - genPt)/(genPt) < -0.05)",\
#                                        "SF = 0.985",                                                                 \
#                                        "L1 p_{T} > 15, Gen p_{T} > 15",                                              \
#                                        "0.985")' 


# root -l -b -q 'makeDistributionPlots.C("SF_0.985_cPt-gt-15_genPt-gt-15_resolution_better_than_-0.05",                \
#                                        "RelVal Electron Gun Pt 2 to 100",                                     \
#                                        "/Users/stephaniekwan/Dropbox/Princeton_G4/Phase2RCT/analyzer/analyzer.root", \
#                                        "(gct_cPt > 15) && (genPt > 15) && ((gct_cPt/0.985 - genPt)/(genPt) > -0.05)",\
#                                        "SF = 0.985",                                                                 \
#                                        "L1 p_{T} > 15, Gen p_{T} > 15",                                              \
#                                        "0.985")'       

# root -l -b -q 'makeDistributionPlots.C("SF_0.985_cPt-gt-15_genPt-gt-15_resolution_between_-0.15_and_-0.05",          \
#                                        "RelVal Electron Gun Pt 2 to 100",                                     \
#                                        "/Users/stephaniekwan/Dropbox/Princeton_G4/Phase2RCT/analyzer/analyzer.root", \
#                                        "(gct_cPt > 15) && (genPt > 15) && ((gct_cPt/0.985 - genPt)/(genPt) < -0.05) && ((gct_cPt/0.985 - genPt)/(genPt) > -0.15)", \
#                                        "SF = 0.985",                                                                 \
#                                        "L1 p_{T} > 15, Gen p_{T} > 15",                                              \
#                                        "0.985")'                        