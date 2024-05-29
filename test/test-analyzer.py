import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process("L1AlgoTest",eras.Phase2C17I13M9)

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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
###113X SinglePion sample 200 PU
#                             fileNames = cms.untracked.vstring('/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/SinglePion_PT0to200/GEN-SIM-DIGI-RAW-MINIAOD/PU200_111X_mcRun4_realistic_T15_v1-v1/120000/02A95838-DEEC-A241-A759-1AECFC6F752F.root'),
###125X SinglePion sample 0 PU
#                            fileNames = cms.untracked.vstring('/store/mc/Phase2Fall22DRMiniAOD/SinglePion_Pt-0To200-gun/GEN-SIM-DIGI-RAW-MINIAOD/noPU_125X_mcRun4_realistic_v2-v1/2550000/14774bf7-4a5b-46a1-b14b-670468f74e77.root'),
###125X SinglePion sample 200 PU
#                            fileNames = cms.untracked.vstring('/store/mc/Phase2Fall22DRMiniAOD/SinglePion_Pt-0To200-gun/GEN-SIM-DIGI-RAW-MINIAOD/PU200_125X_mcRun4_realistic_v2-v1/30000/004d15e3-a12f-4ba9-a2f3-4b7277ffa418.root',
#                            '/store/mc/Phase2Fall22DRMiniAOD/SinglePion_Pt-0To200-gun/GEN-SIM-DIGI-RAW-MINIAOD/PU200_125X_mcRun4_realistic_v2-v1/30000/031844f8-29c9-4911-8cda-0a0a129f01d0.root',
#                            '/store/mc/Phase2Fall22DRMiniAOD/SinglePion_Pt-0To200-gun/GEN-SIM-DIGI-RAW-MINIAOD/PU200_125X_mcRun4_realistic_v2-v1/30000/046bf923-505d-4f28-a226-ef05797b2e83.root',
#                            '/store/mc/Phase2Fall22DRMiniAOD/SinglePion_Pt-0To200-gun/GEN-SIM-DIGI-RAW-MINIAOD/PU200_125X_mcRun4_realistic_v2-v1/30000/0523e9cb-e025-4393-81a0-7f722f2fc90e.root',
#                            '/store/mc/Phase2Fall22DRMiniAOD/SinglePion_Pt-0To200-gun/GEN-SIM-DIGI-RAW-MINIAOD/PU200_125X_mcRun4_realistic_v2-v1/30000/083811e9-2c91-409b-8b1f-4ee2c7ea5934.root',
#                            '/store/mc/Phase2Fall22DRMiniAOD/SinglePion_Pt-0To200-gun/GEN-SIM-DIGI-RAW-MINIAOD/PU200_125X_mcRun4_realistic_v2-v1/30000/0a37df8b-d300-4119-9677-7bd8d098cd4a.root',
#                            '/store/mc/Phase2Fall22DRMiniAOD/SinglePion_Pt-0To200-gun/GEN-SIM-DIGI-RAW-MINIAOD/PU200_125X_mcRun4_realistic_v2-v1/30000/0be55d68-f56a-4578-a4ad-32b25a7f4ecd.root'
#                            ),
###125X VBFH sample 200 PU
#                            fileNames = cms.untracked.vstring('/store/mc/Phase2Fall22DRMiniAOD/VBFHToTauTau_M-125_TuneCP5_14TeV-powheg-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_125X_mcRun4_realistic_v2-v1/30000/00c6c0d5-8072-4710-82bb-1f2ce570ddf4.root',
#                            '/store/mc/Phase2Fall22DRMiniAOD/VBFHToTauTau_M-125_TuneCP5_14TeV-powheg-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_125X_mcRun4_realistic_v2-v1/30000/00e3e318-f1e1-45ae-91c0-82106d2f9d8b.root',
#                            '/store/mc/Phase2Fall22DRMiniAOD/VBFHToTauTau_M-125_TuneCP5_14TeV-powheg-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_125X_mcRun4_realistic_v2-v1/30000/014cc1c6-f0f4-4936-a924-53d8671dbd93.root',
#                            '/store/mc/Phase2Fall22DRMiniAOD/VBFHToTauTau_M-125_TuneCP5_14TeV-powheg-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_125X_mcRun4_realistic_v2-v1/30000/0192404d-3c65-4428-9d82-70240ac70bf6.root',
#                            '/store/mc/Phase2Fall22DRMiniAOD/VBFHToTauTau_M-125_TuneCP5_14TeV-powheg-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_125X_mcRun4_realistic_v2-v1/30000/02fbb2cd-be62-407c-bf92-79491f90f0b2.root',
#                            '/store/mc/Phase2Fall22DRMiniAOD/VBFHToTauTau_M-125_TuneCP5_14TeV-powheg-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_125X_mcRun4_realistic_v2-v1/30000/04045971-9762-4d0c-aef9-d6076a844e14.root',
#                            '/store/mc/Phase2Fall22DRMiniAOD/VBFHToTauTau_M-125_TuneCP5_14TeV-powheg-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_125X_mcRun4_realistic_v2-v1/30000/04f383b5-c2e9-4a6a-b417-7dec76b5987d.root',
#                            '/store/mc/Phase2Fall22DRMiniAOD/VBFHToTauTau_M-125_TuneCP5_14TeV-powheg-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_125X_mcRun4_realistic_v2-v1/30000/06fcf2e5-3d62-4588-8b53-a36fb000c645.root',
#                            ),
###131X VBFH sample 0 PU
#                            fileNames = cms.untracked.vstring('/store/mc/Phase2Spring23DIGIRECOMiniAOD/VBFHToTauTau_M-125_TuneCP5_14TeV-powheg-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/noPU_131X_mcRun4_realistic_v5-v1/2530000/006f6626-cf25-4030-9835-cd507d7fcbec.root'),
###131X VBFH sample 200 PU
                            fileNames = cms.untracked.vstring('/store/mc/Phase2Spring23DIGIRECOMiniAOD/VBFHToTauTau_M-125_TuneCP5_14TeV-powheg-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_Trk1GeV_131X_mcRun4_realistic_v5-v1/30000/005aa288-45c5-4d82-b0d8-4d5849235810.root',
                            '/store/mc/Phase2Spring23DIGIRECOMiniAOD/VBFHToTauTau_M-125_TuneCP5_14TeV-powheg-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_Trk1GeV_131X_mcRun4_realistic_v5-v1/30000/00a27704-5eb4-4b20-81de-8af0d9f09645.root',
                            '/store/mc/Phase2Spring23DIGIRECOMiniAOD/VBFHToTauTau_M-125_TuneCP5_14TeV-powheg-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_Trk1GeV_131X_mcRun4_realistic_v5-v1/30000/071b9cb9-5dd8-4b49-9c47-6e49def34c72.root',),
###Victor
                            #fileNames = cms.untracked.vstring('/store/mc/Phase2Fall22DRMiniAOD/QCD_Pt-15To3000_TuneCP5_Flat_14TeV-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_125X_mcRun4_realistic_v2-v1/40000/00245bdc-37fc-4361-9eec-2f1be8495931.root','/store/mc/Phase2Fall22DRMiniAOD/QCD_Pt-15To3000_TuneCP5_Flat_14TeV-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_125X_mcRun4_realistic_v2-v1/40000/00696033-2e56-4fca-9836-1bac2b3af414.root', '/store/mc/Phase2Fall22DRMiniAOD/QCD_Pt-15To3000_TuneCP5_Flat_14TeV-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_125X_mcRun4_realistic_v2-v1/40000/01637ff3-d279-4850-820d-04062d0bb591.root', '/store/mc/Phase2Fall22DRMiniAOD/QCD_Pt-15To3000_TuneCP5_Flat_14TeV-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_125X_mcRun4_realistic_v2-v1/40000/01c5da90-ac38-4d8e-856d-590564e49d98.root', '/store/mc/Phase2Fall22DRMiniAOD/QCD_Pt-15To3000_TuneCP5_Flat_14TeV-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_125X_mcRun4_realistic_v2-v1/40000/02862c28-2fbd-4211-a9fd-d0b0183a4d40.root', '/store/mc/Phase2Fall22DRMiniAOD/QCD_Pt-15To3000_TuneCP5_Flat_14TeV-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_125X_mcRun4_realistic_v2-v1/40000/03435f05-aca9-4320-b562-3c7cddab675d.root', '/store/mc/Phase2Fall22DRMiniAOD/QCD_Pt-15To3000_TuneCP5_Flat_14TeV-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_125X_mcRun4_realistic_v2-v1/40000/042b6629-64da-4883-9891-6ed221618999.root', '/store/mc/Phase2Fall22DRMiniAOD/QCD_Pt-15To3000_TuneCP5_Flat_14TeV-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_125X_mcRun4_realistic_v2-v1/40000/04a3d93a-6fe3-4c9c-b6b4-1df7cbbc09f9.root', '/store/mc/Phase2Fall22DRMiniAOD/QCD_Pt-15To3000_TuneCP5_Flat_14TeV-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_125X_mcRun4_realistic_v2-v1/40000/057daaa7-6669-41ee-b92e-b03a8421583d.root', '/store/mc/Phase2Fall22DRMiniAOD/QCD_Pt-15To3000_TuneCP5_Flat_14TeV-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_125X_mcRun4_realistic_v2-v1/40000/05b4eb2f-a130-4c0a-bf9b-7b38dddd6897.root',),
                            #fileNames = cms.untracked.vstring('/store/mc/Phase2Fall22DRMiniAOD/QCD_Pt-15To3000_TuneCP5_Flat_14TeV-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_125X_mcRun4_realistic_v2-v1/40000/00696033-2e56-4fca-9836-1bac2b3af414.root'),
###131X DoubleELectron 200 PU
#                            fileNames = cms.untracked.vstring('/store/mc/Phase2Spring23DIGIRECOMiniAOD/DoubleElectron_FlatPt-1To100-gun/GEN-SIM-DIGI-RAW-MINIAOD/PU200_Trk1GeV_131X_mcRun4_realistic_v5-v1/30000/c56be524-bd6e-4c3c-a858-c8934fe8e295.root',
#                            '/store/mc/Phase2Spring23DIGIRECOMiniAOD/DoubleElectron_FlatPt-1To100-gun/GEN-SIM-DIGI-RAW-MINIAOD/PU200_Trk1GeV_131X_mcRun4_realistic_v5-v1/30000/1907ad07-90bc-4299-b45e-7e4df8f3b3aa.root',
#                            '/store/mc/Phase2Spring23DIGIRECOMiniAOD/DoubleElectron_FlatPt-1To100-gun/GEN-SIM-DIGI-RAW-MINIAOD/PU200_Trk1GeV_131X_mcRun4_realistic_v5-v1/30000/f0e011fa-6ed0-4835-a24c-0ef707a6f56f.root',
#                            '/store/mc/Phase2Spring23DIGIRECOMiniAOD/DoubleElectron_FlatPt-1To100-gun/GEN-SIM-DIGI-RAW-MINIAOD/PU200_Trk1GeV_131X_mcRun4_realistic_v5-v1/30000/3a721e81-1e62-4835-a189-892cdba0a6c1.root',
#                            '/store/mc/Phase2Spring23DIGIRECOMiniAOD/DoubleElectron_FlatPt-1To100-gun/GEN-SIM-DIGI-RAW-MINIAOD/PU200_Trk1GeV_131X_mcRun4_realistic_v5-v1/30000/ba7759ce-b9cc-4102-9302-47f51d1a623a.root',),
###Testing event 66894
#                            fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch///store/mc/Phase2Fall22DRMiniAOD/DoubleElectron_FlatPt-1To100-gun/GEN-SIM-DIGI-RAW-MINIAOD/PU200_125X_mcRun4_realistic_v2-v1/30000/0024ebea-73de-496a-9d75-6f0a7c3b2ba4.root'),
#                            fileNames = cms.untracked.vstring('file:/eos/user/p/pdas/calojet/0024ebea-73de-496a-9d75-6f0a7c3b2ba4.root'),
#                            fileNames = cms.untracked.vstring('/store/mc/Phase2Spring23DIGIRECOMiniAOD/DoubleElectron_FlatPt-1To100-gun/GEN-SIM-DIGI-RAW-MINIAOD/PU200_Trk1GeV_131X_mcRun4_realistic_v5-v1/30000/a2817f6f-93c7-4cff-a391-660b9e24d558.root'),
                            inputCommands = cms.untracked.vstring(
                                "keep *",
                                "drop l1tTkPrimaryVertexs_*_*_*",
                            )
                        )

process.source.eventsToProcess = cms.untracked.VEventRange("1:26906")
#process.source.eventsToProcess = cms.untracked.VEventRange("1:44470","1:44374")

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
process.load('L1Trigger.L1CaloTrigger.l1tPhase2CaloJetEmulator_cfi')
process.load('L1Trigger.L1CaloPhase2Analyzer.l1TCaloEGammaAnalyzer_cfi')

process.pL1EG = cms.Path( process.l1tPhase2L1CaloEGammaEmulator*process.l1tPhase2CaloPFClusterEmulator*process.l1tPhase2CaloJetEmulator*process.l1NtupleProducer )

# output file
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('analyzer.root')
)

#process.Out = cms.OutputModule( "PoolOutputModule",
#    fileName = cms.untracked.string( "phase2L1EGammaAnalyzer_VBFH_113X_SinglePion.root" ),
#    outputCommands = cms.untracked.vstring(
#        "keep *_l1tPhase2L1CaloEGammaEmulator_*_*",
#        "keep *_l1tPhase2CaloPFClusterEmulator_*_*",
#        "keep *_l1tHGCalTowerProducer_*_*",
#    )
#)

#process.end = cms.EndPath( process.Out )

process.schedule = cms.Schedule(process.L1simulation_step, process.pL1EG)

# Multi-threading
process.options.numberOfThreads=cms.untracked.uint32(8)
process.options.numberOfStreams=cms.untracked.uint32(0)
