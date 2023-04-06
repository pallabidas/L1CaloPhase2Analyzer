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
                            fileNames = cms.untracked.vstring("root://eos.grif.fr//eos/grif/cms/grif/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/SinglePion_PT0to200/GEN-SIM-DIGI-RAW-MINIAOD/NoPU_111X_mcRun4_realistic_T15_v1-v1/120000/6930E519-DD17-DD40-BE77-41A28DEF9278.root",
                                                              "root://eos.grif.fr//eos/grif/cms/grif/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/SinglePion_PT0to200/GEN-SIM-DIGI-RAW-MINIAOD/NoPU_111X_mcRun4_realistic_T15_v1-v1/120000/A8A57661-B186-C144-AE76-F02CF25EAF61.root",
                                                              "root://eos.grif.fr//eos/grif/cms/grif/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/SinglePion_PT0to200/GEN-SIM-DIGI-RAW-MINIAOD/NoPU_111X_mcRun4_realistic_T15_v1-v1/130000/26B91F7F-C399-074E-9097-B9D87869578D.root",
                                                              "root://eos.grif.fr//eos/grif/cms/grif/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/SinglePion_PT0to200/GEN-SIM-DIGI-RAW-MINIAOD/NoPU_111X_mcRun4_realistic_T15_v1-v1/130000/B4AB0604-1D7F-E549-B5BF-8D90ED66AEC7.root",
                                                              "root://eos.grif.fr//eos/grif/cms/grif/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/SinglePion_PT0to200/GEN-SIM-DIGI-RAW-MINIAOD/NoPU_111X_mcRun4_realistic_T15_v1-v1/270000/75C7ACCC-986B-6440-A933-E1EB0B5494F3.root",
                                                              "root://eos.grif.fr//eos/grif/cms/grif/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/SinglePion_PT0to200/GEN-SIM-DIGI-RAW-MINIAOD/NoPU_111X_mcRun4_realistic_T15_v1-v1/280000/0E148955-73EE-D64E-8A13-DB76DC7159E7.root",
                                                              "root://eos.grif.fr//eos/grif/cms/grif/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/SinglePion_PT0to200/GEN-SIM-DIGI-RAW-MINIAOD/NoPU_111X_mcRun4_realistic_T15_v1-v1/280000/2A52111A-789F-5E46-83AF-01ED60229539.root"
                                                              # root://cmsxrootd.fnal.gov///store/relval/CMSSW_10_6_0_patch2/RelValElectronGunPt2To100/GEN-SIM-DIGI-RAW/106X_upgrade2023_realistic_v3_2023D41noPU-v1/10000/190EDE9F-770B-174A-8BA6-F7814FC67FD4.root,
                                                              # root://cmsxrootd.fnal.gov///store/relval/CMSSW_10_6_0_patch2/RelValElectronGunPt2To100/GEN-SIM-DIGI-RAW/106X_upgrade2023_realistic_v3_2023D41noPU-v1/10000/283255C6-1E20-6F48-8B8B-31E6A62BD48D.root,
                                                              # root://cmsxrootd.fnal.gov///store/relval/CMSSW_10_6_0_patch2/RelValElectronGunPt2To100/GEN-SIM-DIGI-RAW/106X_upgrade2023_realistic_v3_2023D41noPU-v1/10000/71C02E39-ED72-054B-871F-6B1FD1A1C14A.root,
                                                              # root://cmsxrootd.fnal.gov///store/relval/CMSSW_10_6_0_patch2/RelValElectronGunPt2To100/GEN-SIM-DIGI-RAW/106X_upgrade2023_realistic_v3_2023D41noPU-v1/10000/8B75BCAF-FF0C-094C-AB40-08F104148BC0.root, 
                                                              # root://cmsxrootd.fnal.gov///store/relval/CMSSW_10_6_0_patch2/RelValElectronGunPt2To100/GEN-SIM-DIGI-RAW/106X_upgrade2023_realistic_v3_2023D41noPU-v1/10000/A22AA5DB-ACEE-1140-B081-104F634079A1.root 
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

process.load('L1Trigger.L1CaloTrigger.l1tPhase2L1CaloEGammaEmulator_cfi')
process.load('L1Trigger.L1CaloTrigger.l1tPhase2CaloPFClusterEmulator_cfi')
process.load('L1Trigger.L1CaloPhase2Analyzer.l1TCaloEGammaAnalyzer_cfi')

process.pL1EG = cms.Path( process.l1tPhase2L1CaloEGammaEmulator*process.l1tPhase2CaloPFClusterEmulator*process.l1NtupleProducer )

# output file
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('analyzer_singlepion.root')
)

process.Out = cms.OutputModule( "PoolOutputModule",
    fileName = cms.untracked.string( "phase2L1EGammaAnalyzer_singlepion.root" ),
    outputCommands = cms.untracked.vstring(
        "keep *_l1tPhase2L1CaloEGammaEmulator_*_*",
        "keep *_l1tPhase2CaloPFClusterEmulator_*_*",
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
