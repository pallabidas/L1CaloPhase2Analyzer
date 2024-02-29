/*
 *  \file L1TCaloEGammaAnalyzer.cc
 *  Authors S. Kwan, P. Das, I. Ojalvo
 */

// system include files
#include <ap_int.h>
#include <array>
#include <cmath>
// #include <cstdint>
#include <iostream>
#include <fstream>
#include <memory>

// user include files
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CalibFormats/CaloTPG/interface/CaloTPGTranscoder.h"
#include "CalibFormats/CaloTPG/interface/CaloTPGRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
#include "Geometry/HcalTowerAlgo/interface/HcalTrigTowerGeometry.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/L1THGCal/interface/HGCalTower.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"


// ECAL TPs
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"

// HCAL TPs
#include "DataFormats/HcalDigi/interface/HcalTriggerPrimitiveDigi.h"

// Output tower collection
#include "DataFormats/L1TCalorimeterPhase2/interface/CaloCrystalCluster.h"
#include "DataFormats/L1TCalorimeterPhase2/interface/CaloTower.h"
#include "DataFormats/L1TCalorimeterPhase2/interface/CaloPFCluster.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"

#include "L1Trigger/L1CaloTrigger/interface/ParametricCalibration.h"
#include "L1Trigger/L1TCalorimeter/interface/CaloTools.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "L1Trigger/L1CaloPhase2Analyzer/interface/L1TCaloEGammaAnalyzer.h"
#include "DataFormats/Math/interface/deltaR.h"


// ECAL propagation
#include "CommonTools/BaseParticlePropagator/interface/BaseParticlePropagator.h"
#include "CommonTools/BaseParticlePropagator/interface/RawParticle.h"

float etaValues[95] = {-5.2665, -5.1155, -4.92125, -4.71475, -4.53875, -4.36375, -4.1895, -4.014, -3.83875, -3.664, -3.489, -3.314, -3.045, -2.958, -2.871, -2.784, -2.697, -2.61, -2.523, -2.436, -2.349, -2.262, -2.175, -2.088, -2.001, -1.914, -1.827, -1.74, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218, -1.131, -1.044, -0.957, -0.87, -0.783, -0.696, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, -0.087, 0, 0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.87, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479, 1.566, 1.653, 1.74, 1.827, 1.914, 2.001, 2.088, 2.175, 2.262, 2.349, 2.436, 2.523, 2.61, 2.697, 2.784, 2.871, 2.958, 3.045, 3.314, 3.489, 3.664, 3.83875, 4.014, 4.1895, 4.36375, 4.53875, 4.71475, 4.92125, 5.1155, 5.2665};

float phiValues[73] =
    {-3.142, -3.054, -2.967, -2.880, -2.793, -2.705, -2.618, -2.531, -2.443, -2.356, -2.269, -2.182, -2.094, -2.007, -1.920, -1.833, -1.745, -1.658, -1.571, -1.484, -1.396, -1.309, -1.222, -1.134, -1.047, -0.960, -0.873, -0.785, -0.698, -0.611, -0.524, -0.436, -0.349, -0.262, -0.175, -0.087,
     0.000, 0.087, 0.175, 0.262, 0.349, 0.436, 0.524, 0.611, 0.698, 0.785, 0.873, 0.960, 1.047, 1.134, 1.222, 1.309, 1.396, 1.484, 1.571, 1.658, 1.745, 1.833, 1.920, 2.007, 2.094, 2.182, 2.269, 2.356, 2.443, 2.531, 2.618, 2.705, 2.793, 2.880, 2.967, 3.054, 3.142};

using namespace edm;
//using std::cout;
//using std::endl;
//using std::vector;

L1TCaloEGammaAnalyzer::L1TCaloEGammaAnalyzer( const ParameterSet & cfg ) :
  decoderToken_(esConsumes<CaloTPGTranscoder, CaloTPGRecord>(edm::ESInputTag("", ""))),
  caloGeometryToken_(esConsumes<CaloGeometry, CaloGeometryRecord>(edm::ESInputTag("", ""))),
  hbTopologyToken_(esConsumes<HcalTopology, HcalRecNumberingRecord>(edm::ESInputTag("", ""))),
  ecalSrc_(consumes<EcalEBTrigPrimDigiCollection>(cfg.getParameter<edm::InputTag>("ecalDigis"))),
  hcalSrc_(consumes<HcalTrigPrimDigiCollection>(cfg.getParameter<edm::InputTag>("hcalDigis"))),
  rctClustersSrc_(consumes<l1tp2::CaloCrystalClusterCollection>(cfg.getParameter<edm::InputTag>("rctClusters"))),
  gctClustersSrc_(consumes<l1tp2::CaloCrystalClusterCollection>(cfg.getParameter<edm::InputTag>("gctClusters"))),
  rctTowersSrc_(consumes<l1tp2::CaloTowerCollection>(cfg.getParameter<edm::InputTag>("rctTowers"))),
  gctTowersSrc_(consumes<l1tp2::CaloTowerCollection>(cfg.getParameter<edm::InputTag>("gctTowers"))),
  caloPFClustersSrc_(consumes<l1tp2::CaloPFClusterCollection>(cfg.getParameter<edm::InputTag>("PFclusters"))),
  hgcalTowersSrc_(consumes<l1t::HGCalTowerBxCollection>(cfg.getParameter<edm::InputTag>("L1HgcalTowersInputTag"))),
  hfTowersSrc_(consumes<HcalTrigPrimDigiCollection>(cfg.getParameter<edm::InputTag>("hcalDigis"))),
  decoderTag_(esConsumes<CaloTPGTranscoder, CaloTPGRecord>(edm::ESInputTag("", ""))),
  caloJetSrc_(consumes<l1tp2::Phase2L1CaloJetCollection>(cfg.getParameter<edm::InputTag>("caloJets"))),
  recoJetSrc_(consumes<vector<pat::Jet>>(cfg.getParameter<edm::InputTag>("recoJets"))),
  genJetSrc_(consumes<vector<reco::GenJet>>(cfg.getParameter<edm::InputTag>("genJets"))),
  genSrc_ (consumes<std::vector<reco::GenParticle> >(cfg.getParameter<edm::InputTag>("genParticles")))
{
    folderName_          = cfg.getUntrackedParameter<std::string>("folderName");
    h2L1Towers = tfs_->make<TH2F>("h2L1Towers", "Event Display", 94, etaValues, 72, phiValues);
    h2HgcalTowers = tfs_->make<TH2F>("h2HgcalTowers", "Event Display", 94, etaValues, 72, phiValues);
    calo_jet_pt = tfs_->make<TH1F>("calo_jet_pt", "calo_jet_pt", 40, 0., 400.);
    calo_jet_eta = tfs_->make<TH1F>("calo_jet_eta", "calo_jet_eta", 40, -5., 5.);
    calo_jet_phi = tfs_->make<TH1F>("calo_jet_phi", "calo_jet_phi", 40, -M_PI, M_PI);
    reco_jet_pt = tfs_->make<TH1F>("reco_jet_pt", "reco_jet_pt", 40, 0., 400.);
    reco_jet_eta = tfs_->make<TH1F>("reco_jet_eta", "reco_jet_eta", 40, -5., 5.);
    reco_jet_phi = tfs_->make<TH1F>("reco_jet_phi", "reco_jet_phi", 40, -M_PI, M_PI);

    displayTree = tfs_->make<TTree>("displayTree", "Event Display Tree");

    displayTree->Branch("run",    &run,     "run/I");
    displayTree->Branch("lumi",   &lumi,    "lumi/I");
    displayTree->Branch("event",  &event,   "event/I");
    displayTree->Branch("nvtx",   &nvtx,    "nvtx/I");
    
    ////putting bufsize at 32000 and changing split level to 0 so that the branch isn't split into multiple branches
    displayTree->Branch("rctClusters", "vector<TLorentzVector>", &rctClusters, 32000, 0); 
    displayTree->Branch("rctTowers",   "vector<TLorentzVector>", &rctTowers, 32000, 0);
    displayTree->Branch("hcalTPGs", "vector<TLorentzVector>", &allHcalTPGs, 32000, 0); 
    displayTree->Branch("ecalTPGs", "vector<TLorentzVector>", &allEcalTPGs, 32000, 0); 
    displayTree->Branch("hgcalTowers", "vector<TLorentzVector>", &allHgcalTowers, 32000, 0);
    displayTree->Branch("hfTowers", "vector<TLorentzVector>", &allHfTowers, 32000, 0);
    //displayTree->Branch("hgcal_ieta", "vector<int>", &hgcal_ieta, 32000, 0);
    //displayTree->Branch("hgcal_iphi", "vector<int>", &hgcal_iphi, 32000, 0);

    displayTree->Branch("gctTowers",   "vector<TLorentzVector>", &gctTowers, 32000, 0);
    displayTree->Branch("caloPFClusters", "vector<TLorentzVector>", &caloPFClusters, 32000, 0);
    displayTree->Branch("offlineJets", "vector<TLorentzVector>", &offlineJets, 32000, 0);
    displayTree->Branch("genJets", "vector<TLorentzVector>", &genJets, 32000, 0);
    displayTree->Branch("gctCaloJets", "vector<TLorentzVector>", &gctCaloJets, 32000, 0);

    efficiencyTree = tfs_->make<TTree>("efficiencyTree", "Efficiency Tree");
    
    efficiencyTree->Branch("run",    &run,     "run/I");
    efficiencyTree->Branch("lumi",   &lumi,    "lumi/I");
    efficiencyTree->Branch("event",  &event,   "event/I");
    efficiencyTree->Branch("nvtx",   &nvtx,         "nvtx/I");

    // Gen electrons
    efficiencyTree->Branch("genPt",  &genPt,  "genPt/D");
    efficiencyTree->Branch("genEta", &genEta, "genEta/D");
    efficiencyTree->Branch("genPhi", &genPhi, "genPhi/D");

    // The RCT emulator cluster that was matched to the gen electron
    efficiencyTree->Branch("rct_cPt",  &rct_cPt,  "rct_cPt/D");
    efficiencyTree->Branch("rct_cEta", &rct_cEta, "rct_cEta/D");
    efficiencyTree->Branch("rct_cPhi", &rct_cPhi, "rct_cPhi/D");
    efficiencyTree->Branch("rct_deltaR", &rct_deltaR, "rct_deltaR/D");
    efficiencyTree->Branch("rct_et2x5", &rct_et2x5, "rct_et2x5/D");
    efficiencyTree->Branch("rct_et5x5", &rct_et5x5, "rct_et5x5/D");

    // The GCT emulator cluster that was matched to the gen electron
    efficiencyTree->Branch("gct_cPt",  &gct_cPt,  "gct_cPt/D");
    efficiencyTree->Branch("gct_cEta", &gct_cEta, "gct_cEta/D");
    efficiencyTree->Branch("gct_cPhi", &gct_cPhi, "gct_cPhi/D");
    efficiencyTree->Branch("gct_deltaR", &gct_deltaR, "gct_deltaR/D");
    efficiencyTree->Branch("gct_et2x5", &gct_et2x5, "gct_et2x5/D");
    efficiencyTree->Branch("gct_et5x5", &gct_et5x5, "gct_et5x5/D");
    efficiencyTree->Branch("gct_iso",   &gct_iso,   "gct_iso/D");
    efficiencyTree->Branch("gct_is_ss", &gct_is_ss, "gct_is_ss/I");
    efficiencyTree->Branch("gct_is_looseTkss", &gct_is_looseTkss, "gct_is_looseTkss/I");
    efficiencyTree->Branch("gct_is_iso", &gct_is_iso, "gct_is_iso/I");
    efficiencyTree->Branch("gct_is_looseTkiso", &gct_is_looseTkiso, "gct_is_looseTkiso/I");

    pfEfficiencyTree = tfs_->make<TTree>("pfEfficiencyTree", "Efficiency Tree");
    pfEfficiencyTree->Branch("run",    &run,     "run/I");
    pfEfficiencyTree->Branch("lumi",   &lumi,    "lumi/I");
    pfEfficiencyTree->Branch("event",  &event,   "event/I");
    pfEfficiencyTree->Branch("nvtx",   &nvtx,    "nvtx/I");

    // Gen Pion
    pfEfficiencyTree->Branch("genPionPt",  &genPionPt,   "genPionPt/D");
    pfEfficiencyTree->Branch("genPionEta", &genPionEta,  "genPionEta/D");
    pfEfficiencyTree->Branch("genPionPhi", &genPionPhi,  "genPionPhi/D");
 
    // The 3x3 cluster that was matched to the pion
    pfEfficiencyTree->Branch("pf_cPt",  &pf_cPt,  "pf_cPt/D");
    pfEfficiencyTree->Branch("pf_cEta", &pf_cEta, "pf_cEta/D");
    pfEfficiencyTree->Branch("pf_cPhi", &pf_cPhi, "pf_cPhi/D");
    pfEfficiencyTree->Branch("pf_deltaR", &pf_deltaR, "pf_deltaR/D");

    jetEfficiencyTree = tfs_->make<TTree>("jetEfficiencyTree", "Efficiency Tree");

    jetEfficiencyTree->Branch("run",    &run,     "run/I");
    jetEfficiencyTree->Branch("lumi",   &lumi,    "lumi/I");
    jetEfficiencyTree->Branch("event",  &event,   "event/I");
    jetEfficiencyTree->Branch("nvtx",   &nvtx,    "nvtx/I");

    // Gen jets
    jetEfficiencyTree->Branch("genJetPt",  &genJetPt,   "genJetPt/D");
    jetEfficiencyTree->Branch("genJetEta", &genJetEta,  "genJetEta/D");
    jetEfficiencyTree->Branch("genJetPhi", &genJetPhi,  "genJetPhi/D");

    // The GCT jet that was matched to the gen jet
    jetEfficiencyTree->Branch("gctJet_Pt",  &gctJet_Pt,  "gctJet_Pt/D");
    jetEfficiencyTree->Branch("gctJet_Eta", &gctJet_Eta, "gctJet_Eta/D");
    jetEfficiencyTree->Branch("gctJet_Phi", &gctJet_Phi, "gctJet_Phi/D");
    jetEfficiencyTree->Branch("gctJet_deltaR", &gctJet_deltaR, "gctJet_deltaR/D");

  }

void L1TCaloEGammaAnalyzer::beginJob( const EventSetup & es) {
}

void L1TCaloEGammaAnalyzer::analyze( const Event& evt, const EventSetup& es )
 {

  run = evt.id().run();
  lumi = evt.id().luminosityBlock();
  event = evt.id().event();

  //edm::Handle<l1t::HGCalTowerBxCollection> hgcalTowersHandle;
  edm::Handle<l1tp2::Phase2L1CaloJetCollection> caloJets;
  edm::Handle<vector<pat::Jet>> recoJets;
  edm::Handle<vector<reco::GenJet>> genJetColl;

  edm::Handle<l1tp2::CaloCrystalClusterCollection> rctCaloCrystalClusters;
  edm::Handle<l1tp2::CaloTowerCollection> rctCaloL1Towers;
  
  edm::Handle<l1tp2::CaloCrystalClusterCollection> gctCaloCrystalClusters;
  edm::Handle<l1tp2::CaloTowerCollection> gctCaloL1Towers;

  edm::Handle<l1tp2::CaloPFClusterCollection> PFClusters;
  
  edm::Handle<EcalEBTrigPrimDigiCollection> ecalTPGs;
  edm::Handle<HcalTrigPrimDigiCollection> hcalTPGs;  
  edm::Handle<edm::SortedCollection<HcalTriggerPrimitiveDigi> > hbhecoll;
 
  std::vector<Cluster> rctClustersMatched;
  std::vector<Cluster> gctClustersMatched;
  std::vector<TLorentzVector> pfClustersMatched;
  std::vector<TLorentzVector> gctJetMatched;

  std::map<std::string, float> rctExperimentalParams;
  std::map<std::string, float> gctExperimentalParams;

  rctClusters->clear(); 
  rctClusterInfo->clear();
  rctTowers->clear();
  gctClusters->clear();
  gctClusterInfo->clear();
  gctTowers->clear();
  caloPFClusters->clear();
  offlineJets->clear();
  genJets->clear();
  gctCaloJets->clear();
  allEcalTPGs->clear(); 
  allHcalTPGs->clear();
  allHgcalTowers->clear();
  //hgcal_ieta->clear();
  //hgcal_iphi->clear();

  // Detector geometry
  caloGeometry_ = &es.getData(caloGeometryToken_);
  ebGeometry = caloGeometry_->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);
  hbGeometry = caloGeometry_->getSubdetectorGeometry(DetId::Hcal, HcalBarrel);
  hcTopology_ = &es.getData(hbTopologyToken_);
  HcalTrigTowerGeometry theTrigTowerGeometry(hcTopology_);
  decoder_ = &es.getData(decoderToken_);

  // Barrel GCT towers info
  if(evt.getByToken(gctTowersSrc_, gctCaloL1Towers)){
    for(const auto & gctTower : *gctCaloL1Towers){
      TLorentzVector temp;
      temp.SetPtEtaPhiE(gctTower.ecalTowerEt(), gctTower.towerEta(), gctTower.towerPhi(), gctTower.ecalTowerEt());
      gctTowers->push_back(temp);
      h2L1Towers->Fill(gctTower.towerEta(), gctTower.towerPhi(), gctTower.ecalTowerEt());
    }
  }

  // HGCal info
  edm::Handle<l1t::HGCalTowerBxCollection> hgcalTowersHandle;
  if (!evt.getByToken(hgcalTowersSrc_, hgcalTowersHandle))
    std::cout<<"Failed to get towers from hgcalTowerCollection!"<<std::endl;
  evt.getByToken(hgcalTowersSrc_, hgcalTowersHandle);
  l1t::HGCalTowerBxCollection hgcalTowers;
  hgcalTowers = (*hgcalTowersHandle.product());
  for (auto it = hgcalTowers.begin(0); it != hgcalTowers.end(0); it++) {
    float et = it->etEm()+it->etHad();
    float eta = it->eta();
    float phi = it->phi();
    //int ieta = makeEndcapHwIEta(eta);
    //hgcal_ieta->push_back(ieta);
    //int iphi = makeEndcapHwIPhi(phi);
    //hgcal_iphi->push_back(iphi);
    //std::cout<<ieta<<"\t"<<iphi<<"\t"<<et<<std::endl;
    TLorentzVector temp ;
    temp.SetPtEtaPhiE(et,eta,phi,et);
    allHgcalTowers->push_back(temp);
    h2HgcalTowers->Fill(eta, phi, et);
  }

  // HF info
  const auto& decoder = es.getData(decoderTag_);
  edm::Handle<HcalTrigPrimDigiCollection> hfHandle;
  if (!evt.getByToken(hfTowersSrc_, hfHandle))
    std::cout<<"Failed to get HcalTrigPrimDigi for HF!"<<std::endl;
  evt.getByToken(hfTowersSrc_, hfHandle);
  for (const auto& hit : *hfHandle.product()) {
    if (abs(hit.id().ieta()) < l1t::CaloTools::kHFBegin) continue;
    if (abs(hit.id().ieta()) > l1t::CaloTools::kHFEnd) continue;
    float et = decoder.hcaletValue(hit.id(), hit.t0());
    float eta = l1t::CaloTools::towerEta(hit.id().ieta());
    float phi = l1t::CaloTools::towerPhi(hit.id().ieta(), hit.id().iphi());
    TLorentzVector temp ;
    temp.SetPtEtaPhiE(et,eta,phi,et);
    allHfTowers->push_back(temp);
  }

  if(evt.getByToken(caloJetSrc_, caloJets)){
    for(const auto & caloJet : *caloJets){
      TLorentzVector temp;
      temp.SetPtEtaPhiE(caloJet.jetEt(), caloJet.jetEta(), caloJet.jetPhi(), caloJet.jetEt());
      gctCaloJets->push_back(temp);
      calo_jet_pt->Fill(caloJet.jetEt());
      if(caloJet.jetEt() > 50.){
        calo_jet_eta->Fill(caloJet.jetEta());
        calo_jet_phi->Fill(caloJet.jetPhi());
      }
    }
  }

  if(evt.getByToken(recoJetSrc_, recoJets)){
    for(const auto & recoJet : *recoJets){
      TLorentzVector temp;
      temp.SetPtEtaPhiE(recoJet.pt(), recoJet.eta(), recoJet.phi(), recoJet.et());
      offlineJets->push_back(temp);
      reco_jet_pt->Fill(recoJet.pt());
      if(recoJet.pt() > 50.){
        reco_jet_eta->Fill(recoJet.eta());
        reco_jet_phi->Fill(recoJet.phi());
      }
    }
  }

  if(evt.getByToken(genJetSrc_, genJetColl)){
    for(const auto & genJet : *genJetColl){
      TLorentzVector temp;
      temp.SetPtEtaPhiE(genJet.pt(), genJet.eta(), genJet.phi(), genJet.et());
      genJets->push_back(temp);
    }
  }


  // Get the RCT clusters from the emulator and sort them by pT
  if(evt.getByToken(rctClustersSrc_, rctCaloCrystalClusters)){
    for(const auto & rctCluster : *rctCaloCrystalClusters){

//      std::cout << "RCT Cluster found: pT " << rctCluster.pt()  << ", "
//		<< "eta "                   << rctCluster.eta() << ", "
//		<< "phi "                   << rctCluster.phi() << ", " 
//		<< "et 2x5"                 << rctCluster.e2x5() << ", "
//		<< "et 5x5"                 << rctCluster.e5x5() << ", " 
//                << std::endl;     
 
      Cluster temp ;
      TLorentzVector temp_p4;
      temp_p4.SetPtEtaPhiE(rctCluster.pt(),rctCluster.eta(),rctCluster.phi(),rctCluster.pt());
      temp.p4 = temp_p4;

      temp.et2x5 = rctCluster.e2x5();
      temp.et5x5 = rctCluster.e5x5();

      rctExperimentalParams = rctCluster.getExperimentalParams();
      temp.is_ss         = rctExperimentalParams["standaloneWP_showerShape"];
      temp.is_looseTkss  = rctExperimentalParams["trkMatchWP_showerShape"];

//      std::cout << " with flags: "
//                << "is_ss " << temp.is_ss << ","
//                << "is_looseTkss " << temp.is_looseTkss << ", "

      // Save the 4-vector
      rctClusters->push_back(temp_p4);
      // Save the full cluster info 
      rctClusterInfo->push_back(temp);

    }
  }

  std::sort(rctClusters->begin(), rctClusters->end(), L1TCaloEGammaAnalyzer::comparePt);
  std::sort(rctClusterInfo->begin(), rctClusterInfo->end(), L1TCaloEGammaAnalyzer::compareClusterPt);

  // Get the RCT towers from the emulator 
  if(evt.getByToken(rctTowersSrc_, rctCaloL1Towers)) {
    for (const auto & rctTower : *rctCaloL1Towers){
      TLorentzVector temp;
      float totalEt = rctTower.ecalTowerEt() + rctTower.hcalTowerEt();
      if (totalEt > 0) {
//	 std::cout << "Tower found: ECAL ET " << rctTower.ecalTowerEt()  << ", "
//	 	  << "HCAL ET " << rctTower.hcalTowerEt()  << ", "
//	 	  << "iEta, iPhi " << rctTower.towerIEta() << " " << rctTower.towerIPhi() << ", "
//	 	  << "eta "             << rctTower.towerEta() << ", "
//	 	  << "phi "             << rctTower.towerPhi() << std::endl;
	temp.SetPtEtaPhiE(totalEt,
			  rctTower.towerEta(), rctTower.towerPhi(),
			  totalEt);
	rctTowers->push_back(temp);
      }
    }
  }
  
  // Get the GCT clusters from the emulator, and sort them by pT
  if(evt.getByToken(gctClustersSrc_, gctCaloCrystalClusters)){
    for(const auto & gctCluster : *gctCaloCrystalClusters){
      //fill vector                                                                             
      Cluster temp ;
      TLorentzVector temp_p4;
//      std::cout << "GCT Cluster found: pT " << gctCluster.pt()  << ", "
//                << "eta "               << gctCluster.eta() << ", "
//                << "phi "               << gctCluster.phi() << ", " 
//		<< "iso "               << gctCluster.isolation() << ", " 
//		<< std::endl;
      temp_p4.SetPtEtaPhiE(gctCluster.pt(),gctCluster.eta(),gctCluster.phi(),gctCluster.pt());

      temp.p4 = temp_p4;
      temp.et2x5 = gctCluster.e2x5();  // see https://cmssdt.cern.ch/lxr/source/DataFormats/L1TCalorimeterPhase2/interface/CaloCrystalCluster.h
      temp.et5x5 = gctCluster.e5x5();
      temp.iso   = gctCluster.isolation();

      gctExperimentalParams = gctCluster.getExperimentalParams();
      temp.is_ss         = gctExperimentalParams["standaloneWP_showerShape"];
      temp.is_iso        = gctExperimentalParams["standaloneWP_isolation"];
      temp.is_looseTkss  = gctExperimentalParams["trkMatchWP_showerShape"];
      temp.is_looseTkiso = gctExperimentalParams["trkMatchWP_isolation"];

//      std::cout << " with flags: " 
//		<< "is_ss " << temp.is_ss << ","
//		<< "is_iso " << temp.is_iso << ", "
//		<< "is_looseTkss " << temp.is_looseTkss << ", "
//		<< "is_looseTkiso " << temp.is_looseTkiso << std::endl;
      
      // Save the 4-vector
      gctClusters->push_back(temp_p4);
      // Save the full cluster info
      gctClusterInfo->push_back(temp);
    }
  }
  std::sort(gctClusters->begin(), gctClusters->end(), L1TCaloEGammaAnalyzer::comparePt);
  std::sort(gctClusterInfo->begin(), gctClusterInfo->end(), L1TCaloEGammaAnalyzer::compareClusterPt);

  // Get the PF clusters from the emulator, and sort them by pT
  if(evt.getByToken(caloPFClustersSrc_, PFClusters)) {
    for (const auto & PFCluster : *PFClusters){
//      std::cout << "PF cluster found: ET " << PFCluster.clusterEt()  << ", "
//                << "iEta, iPhi " << PFCluster.clusterIEta() << " " << PFCluster.clusterIPhi() << ", "
//                << "eta, phi " << PFCluster.clusterEta() << " " << PFCluster.clusterPhi() << std::endl;
      TLorentzVector temp;
      temp.SetPtEtaPhiE(PFCluster.clusterEt(), PFCluster.clusterEta(), PFCluster.clusterPhi(), PFCluster.clusterEt());
      caloPFClusters->push_back(temp);
    }
  }

  // get the ECAL inputs (i.e. ECAL crystals)
  if(!evt.getByToken(ecalSrc_, ecalTPGs))
    std::cout<<"ERROR GETTING THE ECAL TPGS"<<std::endl;
  else
    for (const auto& hit : *ecalTPGs.product()) {
      if (hit.encodedEt() > 0)  // hit.encodedEt() returns an int corresponding to 2x the crystal Et
	{
	  // Et is 10 bit, by keeping the ADC saturation Et at 120 GeV it means that you have to divide by 8
	  float et = hit.encodedEt() / 8.;
	  
	  if (et < 0.5)
	    continue;  // keep the 500 MeV ET Cut 
	  
	  auto cell = ebGeometry->getGeometry(hit.id());
	  
	  GlobalVector position=GlobalVector(cell->getPosition().x(), cell->getPosition().y(), cell->getPosition().z());
	  float eta = position.eta();
	  float phi = position.phi();
	  TLorentzVector temp ;
	  temp.SetPtEtaPhiE(et,eta,phi,et); 
	  allEcalTPGs->push_back(temp);
	}
    }
  

  if(!evt.getByToken(hcalSrc_, hcalTPGs))
    std::cout<<"ERROR GETTING THE HCAL TPGS"<<std::endl;
  else
  for (const auto& hit : *hcalTPGs.product()) {
    float et = decoder_->hcaletValue(hit.id(), hit.t0());
    ap_uint<10> encodedEt = hit.t0().compressedEt(); 
    // same thing as SOI_compressedEt() in HcalTriggerPrimitiveDigi.h///
    if (et <= 0)
      continue;
    
    if (!(hcTopology_->validHT(hit.id()))) {
      LogError("Phase2L1CaloEGammaEmulator")
  	<< " -- Hcal hit DetID not present in HCAL Geom: " << hit.id() << std::endl;
      throw cms::Exception("Phase2L1CaloEGammaEmulator");
      continue;
    }
    const std::vector<HcalDetId>& hcId = theTrigTowerGeometry.detIds(hit.id());
    if (hcId.empty()) {
      LogError("Phase2L1CaloEGammaEmulator")
  	<< "Cannot find any HCalDetId corresponding to " << hit.id() << std::endl;
      throw cms::Exception("Phase2L1CaloEGammaEmulator");
      continue;
    }
    if (hcId[0].subdetId() > 1)
      continue;
    GlobalVector hcal_tp_position = GlobalVector(0., 0., 0.);
    for (const auto& hcId_i : hcId) {
      if (hcId_i.subdetId() > 1)
        continue;
      // get the first HCAL TP/ cell
      auto cell = hbGeometry->getGeometry(hcId_i);
      if (cell == nullptr)
  	continue;
      GlobalVector tmpVector = GlobalVector(cell->getPosition().x(), cell->getPosition().y(), cell->getPosition().z());
      hcal_tp_position = tmpVector;
      break;
    }
  
    float eta = hcal_tp_position.eta();
    float phi = hcal_tp_position.phi();
    TLorentzVector temp ;
    temp.SetPtEtaPhiE(et,eta,phi,et);
    allHcalTPGs->push_back(temp);
  }

  // Get genParticles
  edm::Handle<GenParticleCollectionType> genParticleHandle;
  if(!evt.getByToken(genSrc_,genParticleHandle)) std::cout<<"No gen Particles Found "<<std::endl;
  
  std::vector<reco::GenParticle> genElectrons;
  std::vector<reco::GenParticle> genPions;
  
  for (unsigned int i = 0; i< genParticleHandle->size(); i++){
    edm::Ptr<reco::GenParticle> ptr(genParticleHandle, i); 
    if ( (abs(ptr->pdgId()) == 11) && ( abs(ptr->eta()) < 1.4841 )) genElectrons.push_back(*ptr);
    if ( (abs(ptr->pdgId()) == 211) && ( abs(ptr->eta()) < 1.4841 )) genPions.push_back(*ptr);
  }

  //************************************************************************************/
  // ECAL propagation of gen electrons
  //************************************************************************************/
  std::vector<TLorentzVector> propagatedGenElectrons;

  for (auto genElectron : genElectrons) {
    RawParticle particle(genElectron.p4());
    particle.setVertex(genElectron.vertex().x(), genElectron.vertex().y(), genElectron.vertex().z(), 0.);
    if (fabs(genElectron.pdgId())==11) particle.setMass(.511);
    else particle.setMass(0.);
    
    int pdgId = genElectron.pdgId();
    if (pdgId > 0)  particle.setCharge( -1.0 ); 
    if (pdgId < 0)  particle.setCharge( 1.0 ); 

    float field_z = 4;
    BaseParticlePropagator prop(particle, 0., 0., field_z);
    prop.propagateToEcalEntrance();
    if( prop.getSuccess() != 0 ) {
      GlobalPoint ecal_pos(prop.particle().vertex().x(), prop.particle().vertex().y(), prop.particle().vertex().z());
      TLorentzVector corrGenElectron;

      corrGenElectron.SetPtEtaPhiM(prop.particle().Pt(),
				   ecal_pos.eta(),
				   ecal_pos.phi(),
				   prop.particle().mass());
      
      propagatedGenElectrons.push_back(corrGenElectron);
    }
  }

  for (auto genElectron : propagatedGenElectrons) {
    genPt = genElectron.Pt();
    genEta = genElectron.Eta();
    genPhi = genElectron.Phi();

    rctClustersMatched.clear();
    gctClustersMatched.clear();

    rct_cPt    = 0;    rct_cEta   = -999;    rct_cPhi   = -999;
    rct_deltaR = 999;

    // Loop through the RCT clusters which are already sorted in decreasing pT, and check for
    // the first cluster within deltaR < 0.5.  
    for (size_t i = 0; i < rctClusterInfo->size(); ++i) {

      float this_rct_deltaR = reco::deltaR(rctClusterInfo->at(i).p4.Eta(), rctClusterInfo->at(i).p4.Phi(),
					   genElectron.Eta(), genElectron.Phi());

      if (this_rct_deltaR < 0.2) {
	
	TLorentzVector temp_p4;
	temp_p4.SetPtEtaPhiE(rctClusterInfo->at(i).p4.Pt(), rctClusterInfo->at(i).p4.Eta(),
			     rctClusterInfo->at(i).p4.Phi(), rctClusterInfo->at(i).p4.M());
	Cluster temp;
	temp.p4 = temp_p4;
	temp.et2x5 = rctClusterInfo->at(i).et2x5;
	temp.et5x5 = rctClusterInfo->at(i).et5x5;
	temp.is_ss = rctClusterInfo->at(i).is_ss;
	temp.is_looseTkss = rctClusterInfo->at(i).is_looseTkss;
	rctClustersMatched.push_back(temp);
        //if (this_rct_deltaR < rct_deltaR) {  // for closest_l1
        //  rct_deltaR = this_rct_deltaR;
        //  rct_cPt = rctClusterInfo->at(i).p4.Pt();
        //  rct_cEta = rctClusterInfo->at(i).p4.Eta();
        //  rct_cPhi = rctClusterInfo->at(i).p4.Phi();
        //  rct_et2x5 = rctClusterInfo->at(i).et2x5;
        //  rct_et5x5 = rctClusterInfo->at(i).et5x5;
        //}
      }
    }
    
    // For this gen electron, sort the matched clusters by pT, and only save the highest pT one
    std::sort(rctClustersMatched.begin(), rctClustersMatched.end(), L1TCaloEGammaAnalyzer::compareClusterPt);
    if (rctClustersMatched.size() > 0) {
      rct_cPt  = rctClustersMatched.at(0).p4.Pt();
      rct_cEta = rctClustersMatched.at(0).p4.Eta();
      rct_cPhi = rctClustersMatched.at(0).p4.Phi();
      rct_deltaR = reco::deltaR(rct_cEta, rct_cPhi,
				genElectron.Eta(), genElectron.Phi());
      rct_et2x5 = rctClustersMatched.at(0).et2x5;
      rct_et5x5 = rctClustersMatched.at(0).et5x5; 

//      std::cout << "--> Matched RCT cluster " << rct_cPt 
//		<< " eta: " << rct_cEta
//		<< " phi: " << rct_cPhi
//		<< " with genElectron " << genPt
//		<< " eta: " << genEta
//		<< " phi: " << genPhi 
//		<< " et2x5: " << rct_et2x5 
//		<< " et5x5: " << rct_et5x5 << std::endl;
    }
    
    //************************************************************************************/
    // Loop through the GCT cluster 4-vectors, and gen-match. Only take the first 
    // GCT cluster which has deltaR < 0.5.
    //************************************************************************************/ 
    
    gct_cPt   = 0;    gct_cEta   = -999;   gct_cPhi  = -999;
    gct_deltaR = 999;
    gct_iso = 0;
    gct_et2x5 = 0; gct_et5x5 = 0;
    gct_is_ss = 0; gct_is_looseTkss = 0;
    gct_is_iso = 0; gct_is_looseTkiso = 0;

    for (size_t i = 0; i < gctClusterInfo->size(); ++i) {
//      std::cout << " gctClusterInfo pT " << gctClusterInfo->at(i).p4.Pt() 
//		<< " eta "               << gctClusterInfo->at(i).p4.Eta()
//		<< " phi "               << gctClusterInfo->at(i).p4.Phi() << std::endl;
      float this_gct_deltaR = reco::deltaR(gctClusterInfo->at(i).p4.Eta(), gctClusterInfo->at(i).p4.Phi(),
					   genElectron.Eta(), genElectron.Phi());
      if (this_gct_deltaR < 0.2) {
	TLorentzVector temp_p4;
        temp_p4.SetPtEtaPhiE(gctClusterInfo->at(i).p4.Pt(), gctClusterInfo->at(i).p4.Eta(), gctClusterInfo->at(i).p4.Phi(), gctClusterInfo->at(i).p4.M());

	Cluster myTemp;
	myTemp.p4    = temp_p4;
	myTemp.iso   = gctClusterInfo->at(i).iso;
	myTemp.et2x5 = gctClusterInfo->at(i).et2x5;
	myTemp.et5x5 = gctClusterInfo->at(i).et5x5;
	myTemp.is_ss = gctClusterInfo->at(i).is_ss;
	myTemp.is_looseTkss = gctClusterInfo->at(i).is_looseTkss;
	myTemp.is_iso = gctClusterInfo->at(i).is_iso;
	myTemp.is_looseTkiso = gctClusterInfo->at(i).is_looseTkiso;

        gctClustersMatched.push_back(myTemp);
        //if(this_gct_deltaR < gct_deltaR) {  // for closest_l1
        //  gct_deltaR = this_gct_deltaR;
        //  gct_cPt = gctClusterInfo->at(i).p4.Pt();
        //  gct_cEta = gctClusterInfo->at(i).p4.Eta();
        //  gct_cPhi = gctClusterInfo->at(i).p4.Phi();
        //  gct_iso = gctClusterInfo->at(i).iso;
        //  gct_et2x5 = gctClusterInfo->at(i).et2x5;
        //  gct_et5x5 = gctClusterInfo->at(i).et5x5;
        //  gct_is_ss = gctClusterInfo->at(i).is_ss;
        //  gct_is_looseTkss = gctClusterInfo->at(i).is_looseTkss;
        //  gct_is_iso = gctClusterInfo->at(i).is_iso;
        //  gct_is_looseTkiso = gctClusterInfo->at(i).is_looseTkiso;
        //}

      }
    }
	
    // For this gen electron, sort the matched clusters by pT, and only save the highest pT one     
    
    std::sort(gctClustersMatched.begin(), gctClustersMatched.end(), L1TCaloEGammaAnalyzer::compareClusterPt);

    if (gctClustersMatched.size() > 0) {
      gct_cPt  = gctClustersMatched.at(0).p4.Pt();
      gct_cEta = gctClustersMatched.at(0).p4.Eta();
      gct_cPhi = gctClustersMatched.at(0).p4.Phi();
      gct_deltaR = reco::deltaR(gct_cEta, gct_cPhi,
    				genElectron.Eta(), genElectron.Phi());
      gct_iso   = gctClustersMatched.at(0).iso;
      gct_et2x5 = gctClustersMatched.at(0).et2x5;
      gct_et5x5 = gctClustersMatched.at(0).et5x5; 
      gct_is_ss = gctClustersMatched.at(0).is_ss;
      gct_is_looseTkss = gctClustersMatched.at(0).is_looseTkss;
      gct_is_iso = gctClustersMatched.at(0).is_iso;
      gct_is_looseTkiso = gctClustersMatched.at(0).is_looseTkiso;
      
//      std::cout << "--> Matched GCT cluster " << gct_cPt
//    		<< " eta: " << gct_cEta
//    		<< " phi: " << gct_cPhi
//    		<< " with genElectron " << genPt
//    		<< " eta: " << genEta
//    		<< " phi: " << genPhi  << ". "
//		<< " iso: "   << gct_iso
//    		<< " et2x5: " << gct_et2x5 
//		<< " et5x5: " << gct_et5x5 
//		<< " is_ss: " << gct_is_ss 
//		<< " is_looseTkss: " << gct_is_looseTkss 
//		<< " is_iso: " << gct_is_iso 
//		<< " is_looseTkiso: " << gct_is_looseTkiso 
//		<< std::endl;
    }

    //efficiencyTree->Fill();

  } // end of loop over gen electrons

    //************************************************************************************/
    // Loop through the PF cluster 4-vetcors, and gen-match. Only take the first
    // PF cluster which has deltaR < 0.5.
    //************************************************************************************/

  for (auto genPion : genPions) {
    genPionPt = genPion.pt();
    genPionEta = genPion.eta();
    genPionPhi = genPion.phi();
    pf_cPt   = 0;    pf_cEta   = -999;   pf_cPhi  = -999;    pf_deltaR = 999;
    pfClustersMatched.clear();

    for (size_t i = 0; i < caloPFClusters->size(); ++i) {
      float this_pf_deltaR = reco::deltaR(caloPFClusters->at(i).Eta(), caloPFClusters->at(i).Phi(),
                                           genPion.eta(), genPion.phi());
      if (this_pf_deltaR < 0.2) {
        TLorentzVector temp = caloPFClusters->at(i);

        pfClustersMatched.push_back(temp);
        //if (this_pf_deltaR < pf_deltaR) {   // for closest_l1
        //  pf_deltaR = this_pf_deltaR;
        //  pf_cPt = caloPFClusters->at(i).Pt();
        //  pf_cEta = caloPFClusters->at(i).Eta();
        //  pf_cPhi = caloPFClusters->at(i).Phi();
        //}
      }
    }

    // For this gen pion, sort the matched clusters by pT, and only save the highest pT one    

    std::sort(pfClustersMatched.begin(), pfClustersMatched.end(), L1TCaloEGammaAnalyzer::comparePt);

    if (pfClustersMatched.size() > 0) {
      pf_cPt  = pfClustersMatched.at(0).Pt();
      pf_cEta = pfClustersMatched.at(0).Eta();
      pf_cPhi = pfClustersMatched.at(0).Phi();
      pf_deltaR = reco::deltaR(pf_cEta, pf_cPhi,
                                genPion.eta(), genPion.phi());

//      std::cout << "--> Matched PF cluster " << pf_cPt
//                << " eta: " << pf_cEta
//                << " phi: " << pf_cPhi
//                << " with genPion " << genPion.Pt()
//                << " eta: " << genEta
//                << " phi: " << genPhi  << ". "
//                << std::endl;
    }

    //pfEfficiencyTree->Fill();

  } // end of loop over gen pions

  //************************************************************************************/ 
  // Loop through the gen jets and match to the GCT jets
  //************************************************************************************/ 
  //for (auto genJet : genJets) {
  for (size_t j = 0; j < genJets->size(); ++j) {

    genJetPt = genJets->at(j).Pt();
    genJetEta = genJets->at(j).Eta();
    genJetPhi = genJets->at(j).Phi();

    gctJetMatched.clear();

    gctJet_Pt    = 0;    gctJet_Eta   = -999;    gctJet_Phi   = -999;
    gctJet_deltaR = 999;

    for (size_t i = 0; i < gctCaloJets->size(); ++i) {
      float this_jet_deltaR = reco::deltaR(gctCaloJets->at(i).Eta(), gctCaloJets->at(i).Phi(),
                                           genJets->at(j).Eta(), genJets->at(j).Phi());
      if (this_jet_deltaR < 0.4) {
        TLorentzVector temp = gctCaloJets->at(i);

        gctJetMatched.push_back(temp);
      }
    }

    std::sort(gctJetMatched.begin(), gctJetMatched.end(), L1TCaloEGammaAnalyzer::comparePt);

    if (gctJetMatched.size() > 0) {
      gctJet_Pt = gctJetMatched.at(0).Pt();
      gctJet_Eta = gctJetMatched.at(0).Eta();
      gctJet_Phi = gctJetMatched.at(0).Phi();
      gctJet_deltaR = reco::deltaR(gctJet_Eta, gctJet_Phi, genJets->at(j).Eta(), genJets->at(j).Phi());
    }
    //jetEfficiencyTree->Fill();

  } // end of loop over gen jets

  displayTree->Fill();
 
 }




void L1TCaloEGammaAnalyzer::endJob() {
}

L1TCaloEGammaAnalyzer::~L1TCaloEGammaAnalyzer(){
}

DEFINE_FWK_MODULE(L1TCaloEGammaAnalyzer);
