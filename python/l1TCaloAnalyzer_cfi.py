import FWCore.ParameterSet.Config as cms


l1NtupleProducer = cms.EDAnalyzer("L1TCaloAnalyzer",
                                  folderName              = cms.untracked.string("firstFolder"),
                                  genParticles     = cms.InputTag("genParticles", "", "HLT"),
                                  ecalDigis = cms.InputTag("simEcalEBTriggerPrimitiveDigis","","HLT"),
                                  hcalDigis = cms.InputTag("simHcalTriggerPrimitiveDigis","","HLT"),
                                  rctClusters = cms.InputTag("l1tPhase2L1CaloEGammaEmulator", "RCTClusters"),
                                  gctClusters = cms.InputTag("l1tPhase2L1CaloEGammaEmulator", "GCTClusters"),
                                  rctTowers = cms.InputTag("l1tPhase2L1CaloEGammaEmulator", "RCTTowers"),
                                  gctTowers = cms.InputTag("l1tPhase2L1CaloEGammaEmulator", "GCTFullTowers"),
                                  PFclusters              = cms.InputTag("l1tPhase2CaloPFClusterEmulator", "GCTPFCluster"),
                                  L1HgcalTowersInputTag   = cms.InputTag("l1tHGCalTowerProducer","HGCalTowerProcessor",""),
                                  caloJets                = cms.InputTag("l1tPhase2CaloJetEmulator", "GCTJet"),
                                  caloJetsDigis           = cms.InputTag("l1tPhase2CaloJetEmulator", "GCTDigitizedJet"),
                                  egtocorr18 = cms.InputTag("l1tPhase2GCTBarrelToCorrelatorLayer1Emulator", "GCTEmDigiClusters"),
                                  pftocorr18 = cms.InputTag("l1tPhase2GCTBarrelToCorrelatorLayer1Emulator", "GCTHadDigiClusters"),
                                  datatocorr18 = cms.InputTag("l1tPhase2CaloToCorrelatorTMI18","DigitizedCaloToCorrelatorTMI18"),
                                  recoJets                = cms.InputTag("slimmedJets","","RECO"),
                                  genJets                 = cms.InputTag("slimmedGenJets","","RECO")
)
