import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process("L1AlgoTest",eras.Phase2C9)

process.load('Configuration.StandardSequences.Services_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.EventContent.EventContent_cff')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# Dataset: 
#   /RelValElectronGunPt2To100/CMSSW_10_6_0_patch2-106X_upgrade2023_realistic_v3_2023D41noPU-v1/GEN-SIM-DIGI-RAW
# xrdcp root://cmsxrootd.fnal.gov///store/relval/CMSSW_10_6_0_patch2/RelValElectronGunPt2To100/GEN-SIM-DIGI-RAW/106X_upgrade2023_realistic_v3_2023D41noPU-v1/10000/190EDE9F-770B-174A-8BA6-F7814FC67FD4.root RelValElectronGunPt2To100_190EDE9F-770B-174A-8BA6-F7814FC67FD4.root

process.source = cms.Source("PoolSource",
                            #fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/s/skkwan/public/phase2RCT/RelValElectronGunPt2To100_190EDE9F-770B-174A-8BA6-F7814FC67FD4.root'
                            fileNames = cms.untracked.vstring("root://cms-xrd-global.cern.ch://store/mc/Phase2HLTTDRWinter20DIGI/DoubleElectron_FlatPt-1To100/GEN-SIM-DIGI-RAW/NoPU_110X_mcRun4_realistic_v3-v2/250000/0361ACC8-9298-424C-A782-DF403BE238BB.root",
                                                              "root://cms-xrd-global.cern.ch://store/mc/Phase2HLTTDRWinter20DIGI/DoubleElectron_FlatPt-1To100/GEN-SIM-DIGI-RAW/NoPU_110X_mcRun4_realistic_v3-v2/250000/0A9A79BB-665C-9B46-99F2-58CCAD3F1EAA.root",
                                                              "root://cms-xrd-global.cern.ch://store/mc/Phase2HLTTDRWinter20DIGI/DoubleElectron_FlatPt-1To100/GEN-SIM-DIGI-RAW/NoPU_110X_mcRun4_realistic_v3-v2/250000/0C75866A-F558-5F44-BA74-5A1917D292E9.root",
                                                              "root://cms-xrd-global.cern.ch://store/mc/Phase2HLTTDRWinter20DIGI/DoubleElectron_FlatPt-1To100/GEN-SIM-DIGI-RAW/NoPU_110X_mcRun4_realistic_v3-v2/250000/0DE33D33-E96C-EC4A-8473-8634B46E92FE.root"
#                            fileNames = cms.untracked.vstring("root://cms-xrd-global.cern.ch://store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/SinglePion_PT0to200/GEN-SIM-DIGI-RAW-MINIAOD/NoPU_111X_mcRun4_realistic_T15_v1-v1/120000/6930E519-DD17-DD40-BE77-41A28DEF9278.root",
#                                                              "root://cms-xrd-global.cern.ch://store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/SinglePion_PT0to200/GEN-SIM-DIGI-RAW-MINIAOD/NoPU_111X_mcRun4_realistic_T15_v1-v1/120000/A8A57661-B186-C144-AE76-F02CF25EAF61.root",
#                                                              "root://cms-xrd-global.cern.ch://store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/SinglePion_PT0to200/GEN-SIM-DIGI-RAW-MINIAOD/NoPU_111X_mcRun4_realistic_T15_v1-v1/130000/26B91F7F-C399-074E-9097-B9D87869578D.root",
#                                                              "root://cms-xrd-global.cern.ch://store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/SinglePion_PT0to200/GEN-SIM-DIGI-RAW-MINIAOD/NoPU_111X_mcRun4_realistic_T15_v1-v1/130000/B4AB0604-1D7F-E549-B5BF-8D90ED66AEC7.root",
#                                                              "root://cms-xrd-global.cern.ch://store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/SinglePion_PT0to200/GEN-SIM-DIGI-RAW-MINIAOD/NoPU_111X_mcRun4_realistic_T15_v1-v1/270000/75C7ACCC-986B-6440-A933-E1EB0B5494F3.root",
#                                                              "root://cms-xrd-global.cern.ch://store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/SinglePion_PT0to200/GEN-SIM-DIGI-RAW-MINIAOD/NoPU_111X_mcRun4_realistic_T15_v1-v1/280000/0E148955-73EE-D64E-8A13-DB76DC7159E7.root",
#                                                              "root://cms-xrd-global.cern.ch://store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/SinglePion_PT0to200/GEN-SIM-DIGI-RAW-MINIAOD/NoPU_111X_mcRun4_realistic_T15_v1-v1/280000/2A52111A-789F-5E46-83AF-01ED60229539.root"
                                                        ),
                            inputCommands = cms.untracked.vstring(
                                "keep *",
                                "drop l1tTkPrimaryVertexs_*_*_*",
                            )
                        )


# --------------------------------------------------------------------------------------------                                                    
#                                                                                                                                                            
# ----   Run the relevant algorithms

# ---- Global Tag :
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '123X_mcRun4_realistic_v3', '')

# Choose a geometry
process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')

# Add HCAL Transcoder
process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
process.load('CalibCalorimetry.CaloTPG.CaloTPGTranscoder_cfi')


# --------------------------------------------------------------------------------------------
#
# ----    Produce the L1EGCrystal clusters using Emulator

process.load('L1Trigger.L1CaloTrigger.Phase2L1CaloEGammaEmulator_cfi')
process.load('L1Trigger.L1CaloPhase2Analyzer.l1TCaloEGammaAnalyzer_cfi')
process.Phase2L1CaloPFClusterEmulatorProducer = cms.EDProducer("Phase2L1CaloPFClusterEmulator")

process.pL1EG = cms.Path( process.Phase2L1CaloEGammaEmulatorProducer*process.Phase2L1CaloPFClusterEmulatorProducer*process.l1NtupleProducer )

#process.pL1EG = cms.Path( process.Phase2L1CaloEGammaEmulatorProducer*process.l1NtupleProducer )

# output file
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('analyzer_doubleelectron.root')
)

process.Out = cms.OutputModule( "PoolOutputModule",
    fileName = cms.untracked.string( "phase2L1EGammaAnalyzer_doubleelectron.root" ),
    outputCommands = cms.untracked.vstring(
        "keep *_Phase2L1CaloEGammaEmulatorProducer_*_*",
        "keep *_Phase2L1CaloPFClusterEmulatorProducer_*_*",
#        "keep *_TriggerResults_*_*",
#        "keep *_simHcalTriggerPrimitiveDigis_*_*",
#        "keep *_EcalEBTrigPrimProducer_*_*"
    )
)


process.end = cms.EndPath( process.Out )

process.schedule = cms.Schedule(process.pL1EG, process.end)

# Multi-threading
process.options.numberOfThreads=cms.untracked.uint32(8)
process.options.numberOfStreams=cms.untracked.uint32(0)

#dump_file = open("dump_file.py", "w")
#dump_file.write(process.dumpPython())
