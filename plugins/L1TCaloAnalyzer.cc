/*
 *  \file L1TCaloAnalyzer.cc
 *  Authors S. Kwan, P. Das, I. Ojalvo
 */

// system include files
#include <ap_int.h>
#include <array>
#include <cmath>
#include <typeinfo>
// #include <cstdint>
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
//using std::cout;
//using std::endl;
//using std::vector;

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
  egDigitizedToCorrelatorTMI18Src_(consumes<l1tp2::DigitizedClusterCorrelatorCollectionTMI18>(cfg.getParameter<edm::InputTag>("egtocorr18"))),
  pfDigitizedToCorrelatorTMI18Src_(consumes<l1tp2::DigitizedPFClusterCorrelatorCollectionTMI18>(cfg.getParameter<edm::InputTag>("pftocorr18"))),
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
    
    ////putting bufsize at 32000 and changing split level to 0 so that the branch isn't split into multiple branches
    ////vector<TLorentzVector>'s usage is causing a runtime error (no dictionary for stl collection)
    /*
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
    displayTree->Branch("gctCaloJetsDigitized", "vector<TLorentzVector>", &gctCaloJetsDigitized, 32000, 0);
    displayTree->Branch("genTaus",  "vector<TLorentzVector>", &genTaus, 32000, 0);
    displayTree->Branch("genQuarks", "vector<TLorentzVector>", &genQuarks, 32000, 0);
    */
    // displayTree->Branch("gctCaloJetsDigitized_size", &gctCaloJetsDigitized_size, "gctCaloJetsDigitized_size/I");
    // displayTree->Branch("gctCaloJetsDigitized_et",   &gctCaloJetsDigitized_et,   "gctCaloJetsDigitized_et[gctCaloJetsDigitized_size]/F");
    // displayTree->Branch("gctCaloJetsDigitized_eta",  &gctCaloJetsDigitized_eta,  "gctCaloJetsDigitized_eta[gctCaloJetsDigitized_size]/F");
    // displayTree->Branch("gctCaloJetsDigitized_phi",  &gctCaloJetsDigitized_phi,  "gctCaloJetsDigitized_phi[gctCaloJetsDigitized_size]/F");
    
    // displayTree->Branch("gctCaloJetsDigitized_et", "vector<ap_uint<16>>", &gctCaloJetsDigitized_et, 32000, 0);
    // displayTree->Branch("gctCaloJetsDigitized_eta", "vector<ap_int<14>>", &gctCaloJetsDigitized_eta, 32000, 0);
    // displayTree->Branch("gctCaloJetsDigitized_phi", "vector<ap_int<13>>", &gctCaloJetsDigitized_phi, 32000, 0);

    displayTree->Branch("gctCaloJetsDigitized_etFloat", "vector<float>", &gctCaloJetsDigitized_etFloat, 32000, 0);
    displayTree->Branch("gctCaloJetsDigitized_etaFloat", "vector<float>", &gctCaloJetsDigitized_etaFloat, 32000, 0);
    displayTree->Branch("gctCaloJetsDigitized_phiFloat", "vector<float>", &gctCaloJetsDigitized_phiFloat, 32000, 0);

    displayTree->Branch("gctCaloJets_et", "vector<float>", &gctCaloJets_et, 32000, 0);
    displayTree->Branch("gctCaloJets_eta", "vector<float>", &gctCaloJets_eta, 32000, 0);
    displayTree->Branch("gctCaloJets_phi", "vector<float>", &gctCaloJets_phi, 32000, 0);
    
    // displayTree->Branch("gctCaloJets_size", &gctCaloJets_size, "gctCaloJets_size/I");
    // displayTree->Branch("gctCaloJets_et",   &gctCaloJets_et,   "gctCaloJets_et[gctCaloJets_size]/F");
    // displayTree->Branch("gctCaloJets_eta",  &gctCaloJets_eta,  "gctCaloJets_eta[gctCaloJets_size]/F");
    // displayTree->Branch("gctCaloJets_phi",  &gctCaloJets_phi,  "gctCaloJets_phi[gctCaloJets_size]/F");

    

  }

void L1TCaloAnalyzer::beginJob( const EventSetup & es) {
}

void L1TCaloAnalyzer::analyze( const Event& evt, const EventSetup& es )
 {

  run = evt.id().run();
  lumi = evt.id().luminosityBlock();
  event = evt.id().event();

  std::cout << " Analyser " << event << std::endl ;

  //edm::Handle<l1t::HGCalTowerBxCollection> hgcalTowersHandle;
  edm::Handle<l1tp2::Phase2L1CaloJetCollection> caloJets;
  edm::Handle<l1tp2::DigitizedL1CaloJetCollection> caloJetsDigis;
  edm::Handle<l1tp2::DigitizedClusterCorrelatorCollectionTMI18> egtocorr18 ;
  edm::Handle<l1tp2::DigitizedPFClusterCorrelatorCollectionTMI18> pftocorr18 ;
  edm::Handle<l1tp2::DigitizedCaloToCorrelatorCollectionTMI18> datatocorr18 ;

  std::map<std::string, float> rctExperimentalParams;
  std::map<std::string, float> gctExperimentalParams;

  gctCaloJets_et->clear();
  gctCaloJets_eta->clear();
  gctCaloJets_phi->clear();
  // gctCaloJetsDigitized_et->clear();
  // gctCaloJetsDigitized_eta->clear();
  // gctCaloJetsDigitized_phi->clear();
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


  int cntr03pos = 0 ;
  int cntr03neg = 0 ;
  int cntr01pos = 0 ;
  int cntr01neg = 0 ;

  int cntr13pos = 0 ;
  int cntr13neg = 0 ;
  int cntr11pos = 0 ;
  int cntr11neg = 0 ;

  int cntr23pos = 0 ;
  int cntr23neg = 0 ;
  int cntr21pos = 0 ;
  int cntr21neg = 0 ;

  ap_uint<64> mydata = 0 ;
  ap_uint<64> dataToCL1Card0[162] = {0} ;
  ap_uint<64> dataToCL1Card1[162] = {0} ;
  ap_uint<64> dataToCL1Card2[162] = {0} ;

  std::cout << " Information for correlator : " << std::endl ;

  if(evt.getByToken(dataDigitizedToCorrelatorTMI18Src_, datatocorr18)) {
    for(const auto & pf : *datatocorr18){

	  std::cout << " Data Size OK: 1=" << datatocorr18->size() <<  std::endl;

	  const ap_uint<64>* ptr0 = pf.link0();
	  const ap_uint<64>* ptr1 = pf.link1();
	  const ap_uint<64>* ptr2 = pf.link2();

	  //------------ card0
	  //
	  for(int i=0; i<162; i++){
		  if((i>1 && i<17)||(i>40 && i<57)||(i>81 && i<98)||(i>121 && i<138)){
		l1tp2::DigitizedClusterCorrelatorTMI18 cluster = l1tp2::DigitizedClusterCorrelatorTMI18(ptr0[i]) ; 
		if(cluster.pt()>20) {
		std::cout << " Card0 word " << i << " EG pt " << cluster.pt() << " eta " << cluster.eta() << " phi " << cluster.phi() << std::endl ;
			}
		  }
		 else {
		l1tp2::DigitizedPFClusterCorrelatorTMI18 cluster = l1tp2::DigitizedPFClusterCorrelatorTMI18(ptr0[i]) ; 
		if(cluster.pt()>20) {
		std::cout << " Card0 word " << i << " PF pt " << cluster.pt() << " eta " << cluster.eta() << " phi " << cluster.phi() << std::endl ;
			}
		 }
	  }
	  //------------ card1
	  //
	  for(int i=0; i<162; i++){
		  if((i>1 && i<17)||(i>40 && i<57)||(i>81 && i<98)||(i>121 && i<138)){
		l1tp2::DigitizedClusterCorrelatorTMI18 cluster = l1tp2::DigitizedClusterCorrelatorTMI18(ptr1[i]) ; 
		if(cluster.pt()>20) {
		std::cout << " Card1 word " << i << " EG pt " << cluster.pt() << " eta " << cluster.eta() << " phi " << cluster.phi() << std::endl ;
			}
		  }
		 else {
		l1tp2::DigitizedPFClusterCorrelatorTMI18 cluster = l1tp2::DigitizedPFClusterCorrelatorTMI18(ptr1[i]) ; 
		if(cluster.pt()>20) {
		std::cout << " Card1 word " << i << " PF pt " << cluster.pt() << " eta " << cluster.eta() << " phi " << cluster.phi() << std::endl ;
			}
		 }
	  }
	  //----------- card2
	  //
	  for(int i=0; i<162; i++){
		  if((i>1 && i<17)||(i>40 && i<57)||(i>81 && i<98)||(i>121 && i<138)){
		l1tp2::DigitizedClusterCorrelatorTMI18 cluster = l1tp2::DigitizedClusterCorrelatorTMI18(ptr2[i]) ; 
		if(cluster.pt()>20) {
		std::cout << " Card2 word " << i << " EG pt " << cluster.pt() << " eta " << cluster.eta() << " phi " << cluster.phi() << std::endl ;
			}
		  }
		 else {
		l1tp2::DigitizedPFClusterCorrelatorTMI18 cluster = l1tp2::DigitizedPFClusterCorrelatorTMI18(ptr2[i]) ; 
		if(cluster.pt()>20) {
		std::cout << " Card2 word " << i << " PF pt " << cluster.pt() << " eta " << cluster.eta() << " phi " << cluster.phi() << std::endl ;
			}
		 }
	  }
	  //------------
  }}


  if(evt.getByToken(pfDigitizedToCorrelatorTMI18Src_, pftocorr18)) {

//	  std::cout << " PF SizeOK:" << pftocorr18->size() <<  std::endl;

    for(const auto & pf : *pftocorr18){
//       std::cout << "    pt: " << pf.pt() << 
//        "    eta: " << pf.eta() << 
//        "    phi: " << pf.phi() << 
//        "    cardNr: " << pf.cardNumber() << 
//        "    slrNr: " << pf.slrNumber() << 
//        "    etaPos: " << pf.etaPositive() << 
//        "    data " << pf.data() << std::endl;

       mydata = pf.data() ;
//       if(pf.pt() > 20) std::cout << hex << "  pf pt: " << pf.pt() << "  " << pf.data()  << " card " << pf.cardNumber() << " slr " << pf.slrNumber() << " PosEta " << pf.etaPositive() << std::endl ; 

       if(pf.slrNumber() == 3 && pf.etaPositive() == 1) goto slr3posp ;
       if(pf.slrNumber() == 3 && pf.etaPositive() == 0) goto slr3negp ;
       if(pf.slrNumber() == 1 && pf.etaPositive() == 1) goto slr1posp ;
       if(pf.slrNumber() == 1 && pf.etaPositive() == 0) goto slr1negp ;

slr3posp:  ;

       if(pf.cardNumber() == 0 && cntr03pos < 24){
       dataToCL1Card0[17+cntr03pos] = mydata ;  
       cntr03pos++ ;  
	goto fillendp ;
       } 
       if(pf.cardNumber() == 1 && cntr13pos < 24){
       dataToCL1Card1[17+cntr13pos] = mydata ;
       cntr13pos++ ;
	goto fillendp ;
       } 
       if(pf.cardNumber() == 2 && cntr23pos < 24){
       dataToCL1Card2[17+cntr23pos] = mydata ;
       cntr23pos++ ;
	goto fillendp ;
       } 
	goto fillendp ;

slr3negp:  ;

       if(pf.cardNumber() == 0 && cntr03neg < 24){
       dataToCL1Card0[57+cntr03neg] = mydata ;
       cntr03neg++ ;
	goto fillendp ;
       } 
       if(pf.cardNumber() == 1 && cntr13neg < 24){
       dataToCL1Card1[57+cntr13neg] = mydata ;
       cntr13neg++ ;
	goto fillendp ;
       } 
       if(pf.cardNumber() == 2 && cntr23neg < 24){
       dataToCL1Card2[57+cntr23neg] = mydata ;
       cntr23neg++ ;
	goto fillendp ;
       } 
	goto fillendp ;

slr1posp:  ;


       if(pf.cardNumber() == 0 && cntr01pos < 24){
       dataToCL1Card0[81+17+cntr01pos] = mydata ;
       cntr01pos++ ;
	goto fillendp ;
       } 
       if(pf.cardNumber() == 1 && cntr11pos < 24){
       dataToCL1Card1[81+17+cntr11pos] = mydata ;
       cntr11pos++ ;
	goto fillendp ;
       } 
       if(pf.cardNumber() == 2 && cntr21pos < 24){
       dataToCL1Card2[81+17+cntr21pos] = mydata ;
       cntr21pos++ ;
	goto fillendp ;
       } 
	goto fillendp ;

slr1negp:  ;

       if(pf.cardNumber() == 0 && cntr01neg < 24){
       dataToCL1Card0[81+57+cntr01neg] = mydata ;
       cntr01neg++ ;
	goto fillendp ;
       } 
       if(pf.cardNumber() == 1 && cntr11neg < 24){
       dataToCL1Card1[81+57+cntr11neg] = mydata ;
       cntr11neg++ ;
	goto fillendp ;
       } 
       if(pf.cardNumber() == 2 && cntr21neg < 24){
       dataToCL1Card2[81+57+cntr21neg] = mydata ;
       cntr21neg++ ;
	goto fillendp ; 
       } 
       
fillendp:  ;

  }}
  

  cntr03pos = 0 ;
  cntr03neg = 0 ;
  cntr01pos = 0 ;
  cntr01neg = 0 ;

  cntr13pos = 0 ;
  cntr13neg = 0 ;
  cntr11pos = 0 ;
  cntr11neg = 0 ;

  cntr23pos = 0 ;
  cntr23neg = 0 ;
  cntr21pos = 0 ;
  cntr21neg = 0 ;


  if(evt.getByToken(egDigitizedToCorrelatorTMI18Src_, egtocorr18)){
//       std::cout << " EG Size:" << egtocorr18->size() <<  std::endl;
    for(const auto & egDigi : *egtocorr18){
//       std::cout << "    pt: " << egDigi.pt() << 
//        "    eta: " << egDigi.eta() << 
//        "    phi: " << egDigi.phi() << 
//        "    hoe: " << egDigi.hoe() << 
//        "    iso: " << egDigi.iso() << 
//        "    shape: " << egDigi.shape() << 
//        "    wp: " << egDigi.wp() <<
//        "    timing: " << egDigi.timing() << 
//        "    brems: " << egDigi.brems() << 
//        "    cardNr: " << egDigi.cardNumber() << 
//        "    slrNr: " << egDigi.slrNumber() << 
//        "    etaPos: " << egDigi.etaPositive() << 
//        "    data " << egDigi.data() << std::endl;

       mydata = egDigi.data() ;
//       if(egDigi.pt() > 20) std::cout << hex << "  eg pt: " << egDigi.pt() << " card " << egDigi.cardNumber() << " slr " << egDigi.slrNumber() << egDigi.etaPositive() << std::endl ; 

       if(egDigi.slrNumber() == 3 && egDigi.etaPositive() == 1) goto slr3pos ;
       if(egDigi.slrNumber() == 3 && egDigi.etaPositive() == 0) goto slr3neg ;
       if(egDigi.slrNumber() == 1 && egDigi.etaPositive() == 1) goto slr1pos ;
       if(egDigi.slrNumber() == 1 && egDigi.etaPositive() == 0) goto slr1neg ;

slr3pos:  ;

       if(egDigi.cardNumber() == 0 && cntr03pos < 16){
       dataToCL1Card0[1+cntr03pos] = mydata ;
       cntr03pos++ ;
	goto fillend ;
       } 
       if(egDigi.cardNumber() == 1 && cntr13pos < 16){
       dataToCL1Card1[1+cntr13pos] = mydata ;
       cntr13pos++ ;
	goto fillend ;
       } 
       if(egDigi.cardNumber() == 2 && cntr23pos < 16){
       dataToCL1Card2[1+cntr23pos] = mydata ;
       cntr23pos++ ;
	goto fillend ;
       } 

slr3neg:  ;

       if(egDigi.cardNumber() == 0 && cntr03neg < 16){
       dataToCL1Card0[41+cntr03neg] = mydata ;
       cntr03neg++ ;
	goto fillend ;
       } 
       if(egDigi.cardNumber() == 1 && cntr13neg < 16){
       dataToCL1Card1[41+cntr13neg] = mydata ;
       cntr13neg++ ;
	goto fillend ;
       } 
       if(egDigi.cardNumber() == 2 && cntr23neg < 16){
       dataToCL1Card2[41+cntr23neg] = mydata ;
       cntr23neg++ ;
	goto fillend ;
       } 

slr1pos:  ;


       if(egDigi.cardNumber() == 0 && cntr01pos < 16){
       dataToCL1Card0[81+1+cntr01pos] = mydata ;
       cntr01pos++ ;
	goto fillend ;
       } 
       if(egDigi.cardNumber() == 1 && cntr11pos < 16){
       dataToCL1Card1[81+1+cntr11pos] = mydata ;
       cntr11pos++ ;
	goto fillend ;
       } 
       if(egDigi.cardNumber() == 2 && cntr21pos < 16){
       dataToCL1Card2[81+1+cntr21pos] = mydata ;
       cntr21pos++ ;
	goto fillend ;
       } 

slr1neg:  ;

       if(egDigi.cardNumber() == 0 && cntr01neg < 16){
       dataToCL1Card0[81+41+cntr01neg] = mydata ;
       cntr01neg++ ;
	goto fillend ;
       } 
       if(egDigi.cardNumber() == 1 && cntr11neg < 16){
       dataToCL1Card1[81+41+cntr11neg] = mydata ;
       cntr11neg++ ;
	goto fillend ;
       } 
       if(egDigi.cardNumber() == 2 && cntr21neg < 16){
       dataToCL1Card2[81+41+cntr21neg] = mydata ;
       cntr21neg++ ;
	goto fillend ;
       } 
fillend:  ;
  }}


  displayTree->Fill();
 
 }




void L1TCaloAnalyzer::endJob() {
}

L1TCaloAnalyzer::~L1TCaloAnalyzer(){
}

DEFINE_FWK_MODULE(L1TCaloAnalyzer);
