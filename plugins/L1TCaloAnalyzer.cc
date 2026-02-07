/*
 *  \file L1TCaloAnalyzer.cc
 *  Authors S. Kwan, P. Das, I. Ojalvo
 */

// system include files
#include <ap_int.h>
#include <array>
#include <cmath>
#include <typeinfo>
#include <iostream>
#include <fstream>
#include <memory>
#include <vector>
#include <TLorentzVector.h>
#ifdef __MAKECINT__
#pragma link C++ class vector<TLorentzVector>+;
#endif

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

#include "L1Trigger/L1CaloPhase2Analyzer/interface/L1TCaloAnalyzer.h"
#include "DataFormats/Math/interface/deltaR.h"


// ECAL propagation
#include "CommonTools/BaseParticlePropagator/interface/BaseParticlePropagator.h"
#include "CommonTools/BaseParticlePropagator/interface/RawParticle.h"

float etaValues[95] = {-5.2665, -5.1155, -4.92125, -4.71475, -4.53875, -4.36375, -4.1895, -4.014, -3.83875, -3.664, -3.489, -3.314, -3.045, -2.958, -2.871, -2.784, -2.697, -2.61, -2.523, -2.436, -2.349, -2.262, -2.175, -2.088, -2.001, -1.914, -1.827, -1.74, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218, -1.131, -1.044, -0.957, -0.87, -0.783, -0.696, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, -0.087, 0, 0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.87, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479, 1.566, 1.653, 1.74, 1.827, 1.914, 2.001, 2.088, 2.175, 2.262, 2.349, 2.436, 2.523, 2.61, 2.697, 2.784, 2.871, 2.958, 3.045, 3.314, 3.489, 3.664, 3.83875, 4.014, 4.1895, 4.36375, 4.53875, 4.71475, 4.92125, 5.1155, 5.2665};

float phiValues[73] =
    {-3.142, -3.054, -2.967, -2.880, -2.793, -2.705, -2.618, -2.531, -2.443, -2.356, -2.269, -2.182, -2.094, -2.007, -1.920, -1.833, -1.745, -1.658, -1.571, -1.484, -1.396, -1.309, -1.222, -1.134, -1.047, -0.960, -0.873, -0.785, -0.698, -0.611, -0.524, -0.436, -0.349, -0.262, -0.175, -0.087,
     0.000, 0.087, 0.175, 0.262, 0.349, 0.436, 0.524, 0.611, 0.698, 0.785, 0.873, 0.960, 1.047, 1.134, 1.222, 1.309, 1.396, 1.484, 1.571, 1.658, 1.745, 1.833, 1.920, 2.007, 2.094, 2.182, 2.269, 2.356, 2.443, 2.531, 2.618, 2.705, 2.793, 2.880, 2.967, 3.054, 3.142};

using namespace edm;

L1TCaloAnalyzer::L1TCaloAnalyzer( const ParameterSet & cfg ) :
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
  caloJetDigitizedSrc_(consumes<l1tp2::DigitizedL1CaloJetCollection>(cfg.getParameter<edm::InputTag>("caloJetsDigis"))),
  egDigitizedToCorrelatorTMI18Src_(consumes<l1tp2::GCTEmDigiClusterCollection>(cfg.getParameter<edm::InputTag>("egtocorr18"))),
  pfDigitizedToCorrelatorTMI18Src_(consumes<l1tp2::GCTHadDigiClusterCollection>(cfg.getParameter<edm::InputTag>("pftocorr18"))),
  dataDigitizedToCorrelatorTMI18Src_(consumes<l1tp2::DigitizedCaloToCorrelatorCollectionTMI18>(cfg.getParameter<edm::InputTag>("datatocorr18"))),
  recoJetSrc_(consumes<vector<pat::Jet>>(cfg.getParameter<edm::InputTag>("recoJets"))),
  genJetSrc_(consumes<vector<reco::GenJet>>(cfg.getParameter<edm::InputTag>("genJets"))),
  genSrc_ (consumes<std::vector<reco::GenParticle> >(cfg.getParameter<edm::InputTag>("genParticles")))
{
    folderName_          = cfg.getUntrackedParameter<std::string>("folderName");

    displayTree = tfs_->make<TTree>("displayTree", "Event Display Tree");

    displayTree->Branch("run",    &run,     "run/I");
    displayTree->Branch("lumi",   &lumi,    "lumi/I");
    displayTree->Branch("event",  &event,   "event/I");
    
  }

void L1TCaloAnalyzer::beginJob( const EventSetup & es) {
}

void L1TCaloAnalyzer::analyze( const Event& evt, const EventSetup& es )
 {

  run = evt.id().run();
  lumi = evt.id().luminosityBlock();
  event = evt.id().event();

  std::cout << " Analyser " << event << std::endl ;

  edm::Handle<l1tp2::Phase2L1CaloJetCollection> caloJets;
  edm::Handle<l1tp2::DigitizedL1CaloJetCollection> caloJetsDigis;
  edm::Handle<l1tp2::DigitizedCaloToCorrelatorCollectionTMI18> datatocorr18 ;
  edm::Handle<l1tp2::GCTEmDigiClusterCollection> egtocorr18;
  edm::Handle<l1tp2::GCTHadDigiClusterCollection> pftocorr18;

  std::map<std::string, float> rctExperimentalParams;
  std::map<std::string, float> gctExperimentalParams;

  gctCaloJets_et->clear();
  gctCaloJets_eta->clear();
  gctCaloJets_phi->clear();
  gctCaloJetsDigitized_etFloat->clear();
  gctCaloJetsDigitized_etaFloat->clear();
  gctCaloJetsDigitized_phiFloat->clear();

  if(evt.getByToken(caloJetSrc_, caloJets)){
    // gctCaloJets_size = caloJets->size();
    for(const auto & caloJet : *caloJets){
      // std::cout << "gctCaloJet:" << std::endl;
      // std::cout << "    Et: " << caloJet.jetEt() << std::endl;
      // std::cout << "    Eta: " << caloJet.jetEta() << std::endl;
      // std::cout << "    Phi: " << caloJet.jetPhi() << std::endl;
      gctCaloJets_et->push_back(caloJet.jetEt());
      gctCaloJets_eta->push_back(caloJet.jetEta());
      gctCaloJets_phi->push_back(caloJet.jetPhi());
    }
  }

  
  if(evt.getByToken(caloJetDigitizedSrc_, caloJetsDigis)){
  //   gctCaloJetsDigitized_size = caloJetsDigis->size();
    for(const auto & caloJetDigi : *caloJetsDigis){
   //    std::cout << "gctCaloJetDigitized:" << std::endl;
   //    std::cout << "    Et: " << caloJetDigi.jetEt() << std::endl;
   //    std::cout << "    Eta: " << caloJetDigi.jetEta() << std::endl;
   //    std::cout << "    Phi: " << caloJetDigi.jetPhi() << std::endl;
   //    gctCaloJetsDigitized_et->push_back(caloJetDigi.pt());
   //    gctCaloJetsDigitized_eta->push_back(caloJetDigi.eta());
   //    gctCaloJetsDigitized_phi->push_back(caloJetDigi.phi());
      gctCaloJetsDigitized_etFloat->push_back(caloJetDigi.ptFloat());
      gctCaloJetsDigitized_etaFloat->push_back(caloJetDigi.etaFloat());
      gctCaloJetsDigitized_phiFloat->push_back(caloJetDigi.phiFloat());
    }
  }

  std::cout << " Information for correlator : " << std::endl ;

  if(evt.getByToken(dataDigitizedToCorrelatorTMI18Src_, datatocorr18)) {
    std::cout << " Data Size OK: 3=" << datatocorr18->size() <<  std::endl;

    int iLink = 0;
    for (const auto & pf : *datatocorr18) {

      const l1tp2::GCTDigiClusterLink& linkptr = pf.linkCard();
    
      //------------ card0
      //

      for(int i=0; i<162; i++){
	const auto& varCluster = linkptr[i];
	if((i>1 && i<33) || (i>81 && i<114)){
	  if (auto* em = std::get_if<l1tp2::GCTEmDigiCluster>(&varCluster)) {
	    l1tp2::GCTEmDigiCluster cluster = *em;
	    if(cluster.pt()>20) {
	      std::cout << " Card "<<iLink<< " word " << i << " EG pt " << cluster.pt() << " eta " << cluster.eta() << " phi " << cluster.phi() << std::endl ;
	      if (cluster.clusterRef().isNonnull()) {
                std::cout << "\t ... Access underlying float cluster pT " << cluster.clusterRef()->pt()
                      << " eta, phi " << cluster.clusterRef()->eta() << ", " << cluster.clusterRef()->phi()
                      << std::endl;
	      }
	    }
	  }
	}
	else {
	  if (auto* pf = std::get_if<l1tp2::GCTHadDigiCluster>(&varCluster)) {
	    l1tp2::GCTHadDigiCluster cluster = *pf;
	    if(cluster.pt()>20) {
	      std::cout << " Card "<<iLink<<" word " << i << " PF pt " << cluster.pt() << " eta " << cluster.eta() << " phi " << cluster.phi() << std::endl ;
	      if (cluster.clusterRef().isNonnull()) {
                std::cout << "\t ... Access underlying float cluster pT " << cluster.clusterRef()->clusterEt()
                      << " eta, phi " << cluster.clusterRef()->clusterEta() << ", " << cluster.clusterRef()->clusterPhi()
		      << " ecal ET " << cluster.clusterRef()->ecalEt()
                      << std::endl;
	      }
            }
	  }
	}
      }
      iLink++;
    }
  }
  //------------
   
//  displayTree->Fill();
 
 }




void L1TCaloAnalyzer::endJob() {
}

L1TCaloAnalyzer::~L1TCaloAnalyzer(){
}

DEFINE_FWK_MODULE(L1TCaloAnalyzer);
