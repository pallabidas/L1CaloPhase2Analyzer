import FWCore.ParameterSet.Config as cms


l1NtupleProducer = cms.EDAnalyzer("L1TCaloEGammaAnalyzer",
                                  folderName              = cms.untracked.string("firstFolder"),
                                  genParticles     = cms.InputTag("genParticles", "", "HLT"),
                                  #packedPfCands           = cms.InputTag("packedPFCandidates"),
                                  #pfCands                 = cms.InputTag("particleFlow"),
                                  ecalDigis = cms.InputTag("simEcalEBTriggerPrimitiveDigis","","HLT"),
                                  hcalDigis = cms.InputTag("simHcalTriggerPrimitiveDigis","","HLT"),
                                  rctClusters = cms.InputTag("l1tPhase2L1CaloEGammaEmulator", "RCT"),
                                  gctClusters = cms.InputTag("l1tPhase2L1CaloEGammaEmulator", "GCT"),
                                  PFclusters              = cms.InputTag("l1tPhase2CaloPFClusterEmulator", "GCTPFCluster"),
                                  L1HgcalTowersInputTag   = cms.InputTag("l1tHGCalTowerProducer","HGCalTowerProcessor",""),
                                  caloJets                = cms.InputTag("l1tPhase2CaloJetEmulator", "GCTJet"),
                                  recoJets                = cms.InputTag("slimmedJets","","RECO")
#                                  clusters  = cms.InputTag('L1EGammaClusterEmuProducer')
)
