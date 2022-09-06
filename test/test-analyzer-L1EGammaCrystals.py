import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process("L1AlgoTest",eras.Phase2C4)

process.load('Configuration.StandardSequences.Services_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.EventContent.EventContent_cff')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# Dataset: 
#   /RelValElectronGunPt2To100/CMSSW_10_6_0_patch2-106X_upgrade2023_realistic_v3_2023D41noPU-v1/GEN-SIM-DIGI-RAW
# xrdcp root://cmsxrootd.fnal.gov///store/relval/CMSSW_10_6_0_patch2/RelValElectronGunPt2To100/GEN-SIM-DIGI-RAW/106X_upgrade2023_realistic_v3_2023D41noPU-v1/10000/190EDE9F-770B-174A-8BA6-F7814FC67FD4.root RelValElectronGunPt2To100_190EDE9F-770B-174A-8BA6-F7814FC67FD4.root

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                                #'file:/afs/cern.ch/work/s/skkwan/public/phase2RCT/RelValElectronGunPt2To100_71C02E39-ED72-054B-871F-6B1FD1A1C14A_1_32_3108.root'
                                'file:/afs/cern.ch/work/s/skkwan/public/phase2RCT/RelValElectronGunPt2To100_71C02E39-ED72-054B-871F-6B1FD1A1C14A.root',
                                'file:/afs/cern.ch/work/s/skkwan/public/phase2RCT/RelValElectronGunPt2To100_8B75BCAF-FF0C-094C-AB40-08F104148BC0.root',
                                'file:/afs/cern.ch/work/s/skkwan/public/phase2RCT/RelValElectronGunPt2To100_190EDE9F-770B-174A-8BA6-F7814FC67FD4.root'
# `
                                # 'root://cmsxrootd.fnal.gov///store/relval/CMSSW_10_6_0_patch2/RelValElectronGunPt2To100/GEN-SIM-DIGI-RAW/106X_upgrade2023_realistic_v3_2023D41noPU-v1/10000/190EDE9F-770B-174A-8BA6-F7814FC67FD4.root',
                                                              # 'root://cmsxrootd.fnal.gov///store/relval/CMSSW_10_6_0_patch2/RelValElectronGunPt2To100/GEN-SIM-DIGI-RAW/106X_upgrade2023_realistic_v3_2023D41noPU-v1/10000/283255C6-1E20-6F48-8B8B-31E6A62BD48D.root',
                                                              # 'root://cmsxrootd.fnal.gov///store/relval/CMSSW_10_6_0_patch2/RelValElectronGunPt2To100/GEN-SIM-DIGI-RAW/106X_upgrade2023_realistic_v3_2023D41noPU-v1/10000/71C02E39-ED72-054B-871F-6B1FD1A1C14A.root',
                                                              # 'root://cmsxrootd.fnal.gov///store/relval/CMSSW_10_6_0_patch2/RelValElectronGunPt2To100/GEN-SIM-DIGI-RAW/106X_upgrade2023_realistic_v3_2023D41noPU-v1/10000/8B75BCAF-FF0C-094C-AB40-08F104148BC0.root', 
                                                              # 'root://cmsxrootd.fnal.gov///store/relval/CMSSW_10_6_0_patch2/RelValElectronGunPt2To100/GEN-SIM-DIGI-RAW/106X_upgrade2023_realistic_v3_2023D41noPU-v1/10000/A22AA5DB-ACEE-1140-B081-104F634079A1.root' 
                                                        ),
                            inputCommands = cms.untracked.vstring(
                                "keep *"
                            )
                        )


# --------------------------------------------------------------------------------------------                                                    
#                                                                                                                                                            
# ----   Run the relevant algorithms

# ---- Global Tag :
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '100X_upgrade2023_realistic_v1', '')

# Choose a geometry
process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')

# Add HCAL Transcoder
process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
process.load('CalibCalorimetry.CaloTPG.CaloTPGTranscoder_cfi')


# --------------------------------------------------------------------------------------------
#
# ----    Produce the L1EGCrystal clusters using Emulator

process.load('L1Trigger.L1CaloTrigger.L1EGammaCrystalsEmulatorProducer_cfi')
process.load('L1Trigger.L1CaloPhase2Analyzer.l1EGammaCrystalsProducerAnalyzer_cfi')

process.pL1EG = cms.Path( process.L1EGammaClusterEmuProducer*process.l1NtupleProducer )

# output file
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('analyzer-l1egammaCrystals.root')
)

process.Out = cms.OutputModule( "PoolOutputModule",
    fileName = cms.untracked.string( "l1egammaCrystals.root" ),
    outputCommands = cms.untracked.vstring(
        "keep *_L1EGammaClusterEmuProducer_*_*",
#        "keep *_TriggerResults_*_*",
#        "keep *_simHcalTriggerPrimitiveDigis_*_*",
#        "keep *_EcalEBTrigPrimProducer_*_*"
    )
)


process.end = cms.EndPath( process.Out )

process.schedule = cms.Schedule(process.pL1EG, process.end)

dump_file = open("dump_file.py", "w")
dump_file.write(process.dumpPython())
