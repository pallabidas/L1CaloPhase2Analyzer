import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process("L1AlgoTest",eras.Phase2C9)

process.load('Configuration.StandardSequences.Services_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.EventContent.EventContent_cff')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )


#Search string: `dataset=/TT*/Phase2HLTTDRWinter20DIGI*/GEN-SIM-DIGI-RAW`
# Picking a PU 200 dataset: `dataset=/TT_TuneCP5_14TeV-powheg-pythia8/Phase2HLTTDRWinter20DIGI-PU200_110X_mcRun4_realistic_v3-v2/GEN-SIM-DIGI-RAW`

# One file (10 GB)
#`xrdcp root://cmsxrootd.fnal.gov///store/mc/Phase2HLTTDRWinter20DIGI/TT_TuneCP5_14TeV-powheg-pythia8/GEN-SIM-DIGI-RAW/PU200_110X_mcRun4_realistic_v3-v2/110000/005E74D6-B50E-674E-89E6-EAA9A617B476.root .`

# file dataset=/SinglePion_PT0to200/Phase2HLTTDRWinter20DIGI-NoPU_110X_mcRun4_realistic_v3-v2/GEN-SIM-DIGI-RAW
# xrdcp root://cms-xrd-global.cern.ch///store/mc/Phase2HLTTDRWinter20DIGI/SinglePion_PT0to200/GEN-SIM-DIGI-RAW/NoPU_110X_mcRun4_realistic_v3-v2/50000/0248A46E-D4D8-514D-A6F9-699E32BDF4B6.root SinglePion_PT0to200_PU200_04AC207E-AF58-C04A-9F90-746DDC628248.root

process.source = cms.Source("PoolSource",
#                            fileNames = cms.untracked.vstring('file:/eos/user/s/skkwan/phase2RCTDevel/SinglePion_PT0to200_NoPU_04AC207E-AF58-C04A-9F90-746DDC628248.root'),
#                            fileNames = cms.untracked.vstring('file:/eos/user/s/skkwan/phase2RCTDevel/005E74D6-B50E-674E-89E6-EAA9A617B476.ro ot'),
#                            fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/s/skkwan/public/phase2RCT/TTBar_005E74D6-B50E-674E-89E6-EAA9A617B476.root'),
                            #fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/s/skkwan/public/phase2RCT/RelValElectronGunPt2To100_71C02E39-ED72-054B-871F-6B1FD1A1C14A_1_32_3108.root'),
                            fileNames = cms.untracked.vstring("root://cms-xrd-global.cern.ch://store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/SinglePion_PT0to200/GEN-SIM-DIGI-RAW-MINIAOD/NoPU_111X_mcRun4_realistic_T15_v1-v1/120000/6930E519-DD17-DD40-BE77-41A28DEF9278.root"),
                            inputCommands = cms.untracked.vstring(
                                "keep *",
                                "drop l1tTkPrimaryVertexs_*_*_*",
                            )
                        )
process.source.eventsToProcess = cms.untracked.VEventRange("1:3014")

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
process.load('L1Trigger.L1CaloPhase2Analyzer.l1TEventDisplayGenerator_cfi')

process.pL1EG = cms.Path( process.l1tPhase2L1CaloEGammaEmulator*process.l1tPhase2CaloPFClusterEmulator*process.l1NtupleProducer )

# output file
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('L1EventDisplay_singlepion.root')
)

process.Out = cms.OutputModule( "PoolOutputModule",
    fileName = cms.untracked.string( "phase2L1CaloEGamma_singlepion.root" ),
    outputCommands = cms.untracked.vstring(
        "keep *_l1tPhase2L1CaloEGammaEmulator_*_*",
#        "keep *_TriggerResults_*_*",
#        "keep *_simHcalTriggerPrimitiveDigis_*_*",
#        "keep *_EcalEBTrigPrimProducer_*_*"
    )
)


process.end = cms.EndPath( process.Out )

process.schedule = cms.Schedule(process.pL1EG, process.end)

#dump_file = open("dump_file.py", "w")
#dump_file.write(process.dumpPython())
