import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process("L1AlgoTest",eras.Phase2C9)

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D86Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D86_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring("root://cms-xrd-global.cern.ch://store/mc/Phase2Spring23DIGIRECOMiniAOD/VBFHToTauTau_M-125_TuneCP5_14TeV-powheg-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_Trk1GeV_131X_mcRun4_realistic_v5-v1/30000/005aa288-45c5-4d82-b0d8-4d5849235810.root",
                            "root://cms-xrd-global.cern.ch://store/mc/Phase2Spring23DIGIRECOMiniAOD/VBFHToTauTau_M-125_TuneCP5_14TeV-powheg-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_Trk1GeV_131X_mcRun4_realistic_v5-v1/30000/00a27704-5eb4-4b20-81de-8af0d9f09645.root",
                                                        ),
                            inputCommands = cms.untracked.vstring(
                                "keep *",
                                "drop l1tTkPrimaryVertexs_*_*_*",
                            )
                        )
#process.source.eventsToProcess = cms.untracked.VEventRange("1:290108")

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
#from L1Trigger.L1THGCal.customTowers import custom_towers_energySplit
#process = custom_towers_energySplit(process)

# Run the emulators and analyzer
process.load('L1Trigger.L1CaloTrigger.l1tPhase2L1CaloEGammaEmulator_cfi')
process.load('L1Trigger.L1CaloTrigger.l1tPhase2CaloPFClusterEmulator_cfi')
process.load('L1Trigger.L1CaloTrigger.l1tPhase2CaloJetEmulator_cfi')
process.load('L1Trigger.L1CaloPhase2Analyzer.l1TCaloEGammaAnalyzer_cfi')

process.pL1EG = cms.Path( process.l1tPhase2L1CaloEGammaEmulator*process.l1tPhase2CaloPFClusterEmulator*process.l1tPhase2CaloJetEmulator*process.l1NtupleProducer )

# output file
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('analyzer_VBFH.root')
)

process.Out = cms.OutputModule( "PoolOutputModule",
    fileName = cms.untracked.string( "phase2L1EGammaAnalyzer_VBFH.root" ),
    outputCommands = cms.untracked.vstring(
        "keep *_l1tPhase2L1CaloEGammaEmulator_*_*",
        "keep *_l1tPhase2CaloPFClusterEmulator_*_*",
        "keep *_l1tHGCalTowerProducer_*_*",
    )
)

process.end = cms.EndPath( process.Out )

process.schedule = cms.Schedule(process.L1simulation_step, process.pL1EG)

# Multi-threading
process.options.numberOfThreads=cms.untracked.uint32(8)
process.options.numberOfStreams=cms.untracked.uint32(0)
