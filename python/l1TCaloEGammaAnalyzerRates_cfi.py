import FWCore.ParameterSet.Config as cms


l1NtupleProducer = cms.EDAnalyzer("L1TCaloEGammaAnalyzerRates",
                                  folderName              = cms.untracked.string("firstFolder"),
                                  genParticles     = cms.InputTag("genParticles", "", "HLT"),
                                  #packedPfCands           = cms.InputTag("packedPFCandidates"),
                                  #pfCands                 = cms.InputTag("particleFlow"),
                                  
                                  ecalDigis = cms.InputTag("simEcalEBTriggerPrimitiveDigis","","HLT"),
                                  hcalDigis = cms.InputTag("simHcalTriggerPrimitiveDigis","","HLT"),
                                  rctClusters = cms.InputTag("Phase2L1CaloEGammaEmulatorProducer", "RCT"),
                                  gctClusters = cms.InputTag("Phase2L1CaloEGammaEmulatorProducer", "GCT")
#                                  clusters  = cms.InputTag('L1EGammaClusterEmuProducer')
)
