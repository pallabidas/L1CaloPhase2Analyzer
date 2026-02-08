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
  genSrc_(consumes<std::vector<reco::GenParticle> >(cfg.getParameter<edm::InputTag>("genParticles")))
{
    folderName_          = cfg.getUntrackedParameter<std::string>("folderName");

    displayTree = tfs_->make<TTree>("displayTree", "Event Display Tree");

    displayTree->Branch("run",    &run,     "run/I");
    displayTree->Branch("lumi",   &lumi,    "lumi/I");
    displayTree->Branch("event",  &event,   "event/I");
    ////putting bufsize at 32000 and changing split level to 0 so that the branch isn't split into multiple branches
    displayTree->Branch("hcalTPGs", "vector<TLorentzVector>", &hcalTPGs, 32000, 0); 
    displayTree->Branch("ecalTPGs", "vector<TLorentzVector>", &ecalTPGs, 32000, 0); 
    displayTree->Branch("egClusters",   "vector<TLorentzVector>", &egClusters, 32000, 0);
    displayTree->Branch("pfClusters", "vector<TLorentzVector>", &pfClusters, 32000, 0);
    displayTree->Branch("offlineJets", "vector<TLorentzVector>", &offlineJets, 32000, 0);
    displayTree->Branch("genJets", "vector<TLorentzVector>", &genJets, 32000, 0);
    displayTree->Branch("gctCaloJets", "vector<TLorentzVector>", &gctCaloJets, 32000, 0);
    displayTree->Branch("gctDigiJets", "vector<TLorentzVector>", &gctDigiJets, 32000, 0);
    displayTree->Branch("genEles",  "vector<TLorentzVector>", &genEles, 32000, 0);

    efficiencyTree = tfs_->make<TTree>("efficiencyTree", "Efficiency Tree");

    // Gen electrons
    efficiencyTree->Branch("genPt",  &genPt,  "genPt/D");
    efficiencyTree->Branch("genEta", &genEta, "genEta/D");
    efficiencyTree->Branch("genPhi", &genPhi, "genPhi/D");

    // The EG cluster that was matched to the gen electron
    efficiencyTree->Branch("eg_cPt",  &eg_cPt,  "eg_cPt/D");
    efficiencyTree->Branch("eg_cEta", &eg_cEta, "eg_cEta/D");
    efficiencyTree->Branch("eg_cPhi", &eg_cPhi, "eg_cPhi/D");
    efficiencyTree->Branch("eg_deltaR", &eg_deltaR, "eg_deltaR/D");
    //efficiencyTree->Branch("eg_hoe", &eg_hoe, "eg_hoe/D");
    //efficiencyTree->Branch("eg_shape", &eg_shape, "eg_shape/D");
    //efficiencyTree->Branch("eg_iso",   &eg_iso,   "eg_iso/D");
    //efficiencyTree->Branch("eg_wp", &eg_wp, "eg_wp/I");

    pfEfficiencyTree = tfs_->make<TTree>("pfEfficiencyTree", "Efficiency Tree");

    // Gen jet
    pfEfficiencyTree->Branch("genJetPt",  &genJetPt,   "genJetPt/D");
    pfEfficiencyTree->Branch("genJetEta", &genJetEta,  "genJetEta/D");
    pfEfficiencyTree->Branch("genJetPhi", &genJetPhi,  "genJetPhi/D");
 
    // The PF cluster that was matched to the gen jet
    pfEfficiencyTree->Branch("pf_cPt",  &pf_cPt,  "pf_cPt/D");
    pfEfficiencyTree->Branch("pf_cEta", &pf_cEta, "pf_cEta/D");
    pfEfficiencyTree->Branch("pf_cPhi", &pf_cPhi, "pf_cPhi/D");
    pfEfficiencyTree->Branch("pf_deltaR", &pf_deltaR, "pf_deltaR/D");
    //pfEfficiencyTree->Branch("pf_ecal", &pf_ecal, "pf_ecal/D");

  }

void L1TCaloAnalyzer::beginJob( const EventSetup & es) {
}

void L1TCaloAnalyzer::analyze( const Event& evt, const EventSetup& es ) {

  run = evt.id().run();
  lumi = evt.id().luminosityBlock();
  event = evt.id().event();

  edm::Handle<l1tp2::Phase2L1CaloJetCollection> caloJets;
  edm::Handle<l1tp2::DigitizedL1CaloJetCollection> caloJetsDigis;
  edm::Handle<l1tp2::DigitizedCaloToCorrelatorCollectionTMI18> datatocorr18 ;
  edm::Handle<l1tp2::GCTEmDigiClusterCollection> egtocorr18;
  edm::Handle<l1tp2::GCTHadDigiClusterCollection> pftocorr18;

  edm::Handle<vector<pat::Jet>> recoJets;
  edm::Handle<vector<reco::GenJet>> genJetColl;
  edm::Handle<EcalEBTrigPrimDigiCollection> ecalTPGColl;
  edm::Handle<HcalTrigPrimDigiCollection> hcalTPGColl;
  edm::Handle<edm::SortedCollection<HcalTriggerPrimitiveDigi> > hbhecoll;

  std::map<std::string, float> rctExperimentalParams;
  std::map<std::string, float> gctExperimentalParams;

  hcalTPGs->clear();
  ecalTPGs->clear();
  egClusters->clear();
  pfClusters->clear();
  offlineJets->clear();
  genJets->clear();
  gctCaloJets->clear();
  gctDigiJets->clear();
  genEles->clear();

  // Detector geometry
  caloGeometry_ = &es.getData(caloGeometryToken_);
  ebGeometry = caloGeometry_->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);
  hbGeometry = caloGeometry_->getSubdetectorGeometry(DetId::Hcal, HcalBarrel);
  hcTopology_ = &es.getData(hbTopologyToken_);
  HcalTrigTowerGeometry theTrigTowerGeometry(hcTopology_);
  decoder_ = &es.getData(decoderToken_);

  // get the ECAL inputs (i.e. ECAL crystals)
  if(!evt.getByToken(ecalSrc_, ecalTPGColl))
    std::cout<<"ERROR GETTING THE ECAL TPGS"<<std::endl;
  else
    for (const auto& hit : *ecalTPGColl.product()) {
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
	  ecalTPGs->push_back(temp);
	}
    }
  

  if(!evt.getByToken(hcalSrc_, hcalTPGColl))
    std::cout<<"ERROR GETTING THE HCAL TPGS"<<std::endl;
  else
  for (const auto& hit : *hcalTPGColl.product()) {
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
    hcalTPGs->push_back(temp);
  }

  // Get genParticles
  edm::Handle<GenParticleCollectionType> genParticleHandle;
  if(!evt.getByToken(genSrc_,genParticleHandle)) std::cout<<"No gen Particles Found "<<std::endl;
  
  std::vector<reco::GenParticle> genElectrons;
  
  for (unsigned int i = 0; i< genParticleHandle->size(); i++){
    edm::Ptr<reco::GenParticle> ptr(genParticleHandle, i); 
    if ( (abs(ptr->pdgId()) == 11) && ( abs(ptr->eta()) < 1.4841 )) {
      genElectrons.push_back(*ptr);
      TLorentzVector temp;
      temp.SetPtEtaPhiE(ptr->pt(), ptr->eta(), ptr->phi(), ptr->energy());
      genEles->push_back(temp);
    }
  }

  // ECAL propagation of gen electrons
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

  if(evt.getByToken(caloJetSrc_, caloJets)){
    for(const auto & caloJet : *caloJets){
      TLorentzVector temp;
      temp.SetPtEtaPhiE(caloJet.jetEt(), caloJet.jetEta(), caloJet.jetPhi(), caloJet.jetEt());
      gctCaloJets->push_back(temp);
    }
  }

  if(evt.getByToken(recoJetSrc_, recoJets)){
    for(const auto & recoJet : *recoJets){
      TLorentzVector temp;
      temp.SetPtEtaPhiE(recoJet.pt(), recoJet.eta(), recoJet.phi(), recoJet.et());
      offlineJets->push_back(temp);
    }
  }

  if(evt.getByToken(genJetSrc_, genJetColl)){
    for(const auto & genJet : *genJetColl){
      TLorentzVector temp;
      temp.SetPtEtaPhiE(genJet.pt(), genJet.eta(), genJet.phi(), genJet.et());
      genJets->push_back(temp);
    }
  }

  if(evt.getByToken(caloJetDigitizedSrc_, caloJetsDigis)){
    for(const auto & caloJetDigi : *caloJetsDigis){
      TLorentzVector temp;
      temp.SetPtEtaPhiE(caloJetDigi.ptFloat(), caloJetDigi.etaFloat(), caloJetDigi.phiFloat(), caloJetDigi.ptFloat());
      gctDigiJets->push_back(temp);
    }
  }

  std::cout << " Information for correlator : " << std::endl ;

  if(evt.getByToken(dataDigitizedToCorrelatorTMI18Src_, datatocorr18)) {
    std::cout << " Data Size OK: 3=" << datatocorr18->size() <<  std::endl;

    int iLink = 0;
    for (const auto & pf : *datatocorr18) {

      const l1tp2::GCTDigiClusterLink& linkptr = pf.linkCard();
    
      for(int i=0; i<162; i++){
	const auto& varCluster = linkptr[i];
	if((i>1 && i<33) || (i>81 && i<114)){
	  if (auto* em = std::get_if<l1tp2::GCTEmDigiCluster>(&varCluster)) {
	    l1tp2::GCTEmDigiCluster cluster = *em;
	    if(cluster.pt()>0) {
	      //std::cout << " Card "<<iLink<< " word " << i << " EG pt " << cluster.pt() << " eta " << cluster.eta() << " phi " << cluster.phi() << std::endl ;
	      if (cluster.clusterRef().isNonnull()) {
                //std::cout << "\t ... Access underlying float cluster pT " << cluster.clusterRef()->pt()
                //      << " eta, phi " << cluster.clusterRef()->eta() << ", " << cluster.clusterRef()->phi()
                //      << std::endl;
		TLorentzVector temp;
		temp.SetPtEtaPhiE(cluster.clusterRef()->pt(), cluster.clusterRef()->eta(), cluster.clusterRef()->phi(), cluster.clusterRef()->pt());
		egClusters->push_back(temp);
	      }
	    }
	  }
	}
	else {
	  if (auto* pf = std::get_if<l1tp2::GCTHadDigiCluster>(&varCluster)) {
	    l1tp2::GCTHadDigiCluster cluster = *pf;
	    if(cluster.pt()>0) {
	      //std::cout << " Card "<<iLink<<" word " << i << " PF pt " << cluster.pt() << " eta " << cluster.eta() << " phi " << cluster.phi() << std::endl ;
	      if (cluster.clusterRef().isNonnull()) {
                //std::cout << "\t ... Access underlying float cluster pT " << cluster.clusterRef()->clusterEt()
                //      << " eta, phi " << cluster.clusterRef()->clusterEta() << ", " << cluster.clusterRef()->clusterPhi()
		//      << " ecal ET " << cluster.clusterRef()->ecalEt()
                //      << std::endl;
		TLorentzVector temp;
		temp.SetPtEtaPhiE(cluster.clusterRef()->clusterEt(), cluster.clusterRef()->clusterEta(), cluster.clusterRef()->clusterPhi(), cluster.clusterRef()->clusterEt());
		pfClusters->push_back(temp);
	      }
            }
	  }
	}
      }
      iLink++;
    }
  }
  //------------
   
  displayTree->Fill();

  //------------
  std::sort(egClusters->begin(), egClusters->end(), comparePt);
  std::sort(pfClusters->begin(), pfClusters->end(), comparePt);

  for (auto gen : propagatedGenElectrons) {
    genPt = gen.Pt();
    genEta = gen.Eta();
    genPhi = gen.Phi();
    eg_deltaR = 0.2; eg_cPt = -99.; eg_cEta = -99.; eg_cPhi = -99.;

    for (size_t i = 0; i < egClusters->size(); ++i) {
      float tempDR = reco::deltaR(egClusters->at(i).Eta(), egClusters->at(i).Phi(), genEta, genPhi);
      if (tempDR < eg_deltaR) {
	eg_deltaR = tempDR;
	eg_cPt = egClusters->at(i).Pt();
	eg_cEta = egClusters->at(i).Eta();
	eg_cPhi = egClusters->at(i).Phi();
      }
    }
    efficiencyTree->Fill();
  }

  for (size_t j = 0; j < genJets->size(); ++j) {
    genJetPt = genJets->at(j).Pt();
    genJetEta = genJets->at(j).Eta();
    genJetPhi = genJets->at(j).Phi();
    pf_deltaR = 0.2; pf_cPt = -99.; pf_cEta = -99.; pf_cPhi = -99.;    
    for (size_t i = 0; i < pfClusters->size(); ++i) {
      float tempDR = reco::deltaR(pfClusters->at(i).Eta(), pfClusters->at(i).Phi(), genJetEta, genJetPhi);
      if (tempDR < pf_deltaR) {
        pf_deltaR = tempDR;
        pf_cPt = pfClusters->at(i).Pt();
        pf_cEta = pfClusters->at(i).Eta();
        pf_cPhi = pfClusters->at(i).Phi();
      }
    }
    pfEfficiencyTree->Fill();
  }
 
}

void L1TCaloAnalyzer::endJob() {
}

L1TCaloAnalyzer::~L1TCaloAnalyzer(){
}

DEFINE_FWK_MODULE(L1TCaloAnalyzer);
