import FWCore.ParameterSet.Config as cms


l1NtupleProducer = cms.EDAnalyzer("L1TEventDisplayGenerator",
                                  folderName              = cms.untracked.string("firstFolder"),
                                  #packedPfCands           = cms.InputTag("packedPFCandidates"),
                                  #pfCands                 = cms.InputTag("particleFlow"),
                                  ecalDigis = cms.InputTag("simEcalEBTriggerPrimitiveDigis","","HLT"),
                                  hcalDigis = cms.InputTag("simHcalTriggerPrimitiveDigis","","HLT"),
                                  clusters                = cms.InputTag("Phase2L1CaloEGammaEmulatorProducer", "GCT"),
                                  towers                  = cms.InputTag("Phase2L1CaloEGammaEmulatorProducer", "GCTFullTowers"), 
                                  PFclusters              = cms.InputTag("Phase2L1CaloPFClusterEmulatorProducer", "GCTPFCluster")
#                                  clusters  = cms.InputTag('L1EGammaClusterEmuProducer')
)
