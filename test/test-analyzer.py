import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process("L1AlgoTest",eras.Phase2C17I13M9)

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtendedRun4D110Reco_cff')
process.load('Configuration.Geometry.GeometryExtendedRun4D110_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                                "root://cms-xrd-global.cern.ch://store/mc/Phase2Spring24DIGIRECOMiniAOD/DoublePhoton_FlatPt-1To100-gun/GEN-SIM-DIGI-RAW-MINIAOD/PU200_Trk1GeV_140X_mcRun4_realistic_v4-v2/2560000/4b0a0e07-247b-458e-9540-2473353ac6d8.root",),
#                                "root://cms-xrd-global.cern.ch://store/mc/Phase2Spring24DIGIRECOMiniAOD/DoubleElectron_FlatPt-1To100-gun/GEN-SIM-DIGI-RAW-MINIAOD/PU200_Trk1GeV_140X_mcRun4_realistic_v4-v2/2810000/001ebf5f-b83c-43fc-997f-c2e5ecf1f9dd.root",),
#/store/mc/Phase2Spring24DIGIRECOMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_Trk1GeV_140X_mcRun4_realistic_v4-v1/2810000/fefd9caa-182c-418e-8419-afb40a9a7739.root
#/store/mc/Phase2Spring24DIGIRECOMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_Trk1GeV_140X_mcRun4_realistic_v4-v1/2810000/ff04ba49-7295-492c-b8de-84d61d9836d7.root
#/store/mc/Phase2Spring24DIGIRECOMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_Trk1GeV_140X_mcRun4_realistic_v4-v1/2810000/ff13f279-3532-4511-8a75-2de7bc476540.root
#/store/mc/Phase2Spring24DIGIRECOMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_Trk1GeV_140X_mcRun4_realistic_v4-v1/2810000/ff2b5a9b-f5b5-4c79-8511-40e1d6a1ce66.root
#/store/mc/Phase2Spring24DIGIRECOMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_Trk1GeV_140X_mcRun4_realistic_v4-v1/2810000/ff6b1222-b2c8-47a6-a265-03e5b885b130.root
#/store/mc/Phase2Spring24DIGIRECOMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_Trk1GeV_140X_mcRun4_realistic_v4-v1/2810000/ff8b25af-a043-4c03-ab78-d2d2be95ce5d.root

                            inputCommands = cms.untracked.vstring(
                                "keep *",
                                "drop l1tTkPrimaryVertexs_*_*_*",
                            )
                        )

# Global Tag :
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '131X_mcRun4_realistic_v6', '')

# Add HCAL Transcoder
process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
process.load('CalibCalorimetry.CaloTPG.CaloTPGTranscoder_cfi')

# Run L1 simulation (necessary to get HGCal towers)
process.L1simulation_step = cms.Path(process.SimL1Emulator)

# Add HGCal module energy splitting
from L1Trigger.L1THGCal.customTowers import custom_towers_energySplit
process = custom_towers_energySplit(process)

# Run the emulators and analyzer
process.load('L1Trigger.L1CaloTrigger.l1tPhase2L1CaloEGammaEmulator_cfi')
process.load('L1Trigger.L1CaloTrigger.l1tPhase2CaloPFClusterEmulator_cfi')
process.load('L1Trigger.L1CaloTrigger.l1tPhase2CaloToCorrelatorTMI18_cfi')
process.load('L1Trigger.L1CaloTrigger.l1tPhase2CaloJetEmulator_cfi')
process.load('L1Trigger.L1CaloPhase2Analyzer.l1TCaloAnalyzer_cfi')

process.pL1EG = cms.Path( process.l1tPhase2L1CaloEGammaEmulator*process.l1tPhase2CaloPFClusterEmulator*process.l1tPhase2CaloJetEmulator*process.l1tPhase2CaloToCorrelatorTMI18*process.l1NtupleProducer )

# output file
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('analyzer.root')
)

process.schedule = cms.Schedule(process.L1simulation_step, process.pL1EG)

# Multi-threading
#process.options.numberOfThreads=cms.untracked.uint32(8)
#process.options.numberOfStreams=cms.untracked.uint32(0)
