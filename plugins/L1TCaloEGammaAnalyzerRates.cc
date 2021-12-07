/*
 *  \file L1TCaloEGammaAnalyzerRates.cc
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

// ECAL TPs
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"

// HCAL TPs
#include "DataFormats/HcalDigi/interface/HcalTriggerPrimitiveDigi.h"

// Output tower collection
#include "DataFormats/L1TCalorimeterPhase2/interface/CaloCrystalCluster.h"
#include "DataFormats/L1TCalorimeterPhase2/interface/CaloTower.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"

#include "L1Trigger/L1CaloTrigger/interface/ParametricCalibration.h"
#include "L1Trigger/L1TCalorimeter/interface/CaloTools.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "L1Trigger/L1CaloPhase2Analyzer/interface/L1TCaloEGammaAnalyzerRates.h"
#include "DataFormats/Math/interface/deltaR.h"


// ECAL propagation
#include "CommonTools/BaseParticlePropagator/interface/BaseParticlePropagator.h"
#include "CommonTools/BaseParticlePropagator/interface/RawParticle.h"


using namespace edm;
using std::cout;
using std::endl;
using std::vector;

L1TCaloEGammaAnalyzerRates::L1TCaloEGammaAnalyzerRates( const ParameterSet & cfg ) :
  ecalSrc_(consumes<EcalEBTrigPrimDigiCollection>(cfg.getParameter<edm::InputTag>("ecalDigis"))),
  hcalSrc_(consumes<HcalTrigPrimDigiCollection>(cfg.getParameter<edm::InputTag>("hcalDigis"))),
  rctClustersSrc_(consumes<l1tp2::CaloCrystalClusterCollection >(cfg.getParameter<edm::InputTag>("rctClusters"))),
  gctClustersSrc_(consumes<l1tp2::CaloCrystalClusterCollection >(cfg.getParameter<edm::InputTag>("gctClusters"))),
  rctTowersSrc_(consumes<l1tp2::CaloTowerCollection >(cfg.getParameter<edm::InputTag>("rctClusters"))),
 gctTowersSrc_(consumes<l1tp2::CaloTowerCollection >(cfg.getParameter<edm::InputTag>("gctClusters"))),
  genSrc_ (( cfg.getParameter<edm::InputTag>( "genParticles")))
{
  genToken_ =     consumes<std::vector<reco::GenParticle> >(genSrc_);

    folderName_          = cfg.getUntrackedParameter<std::string>("folderName");
    efficiencyTree = tfs_->make<TTree>("efficiencyTree", "Efficiency Tree");
    
    ////putting bufsize at 32000 and changing split level to 0 so that the branch isn't split into multiple branches
    efficiencyTree->Branch("rctClusters", "vector<TLorentzVector>", &rctClusters, 32000, 0); 
    efficiencyTree->Branch("rctTowers",   "vector<TLorentzVector>", &rctTowers, 32000, 0);
    efficiencyTree->Branch("hcalTPGs", "vector<TLorentzVector>", &allHcalTPGs, 32000, 0); 
    efficiencyTree->Branch("ecalTPGs", "vector<TLorentzVector>", &allEcalTPGs, 32000, 0); 

    efficiencyTree->Branch("gctClusters", "vector<TLorentzVector>", &gctClusters, 32000, 0);
    efficiencyTree->Branch("gctTowers",   "vector<TLorentzVector>", &gctTowers, 32000, 0);
    
    efficiencyTree->Branch("run",    &run,     "run/I");
    efficiencyTree->Branch("lumi",   &lumi,    "lumi/I");
    efficiencyTree->Branch("event",  &event,   "event/I");

    // Gen electron (only highest pT one per event)
    efficiencyTree->Branch("genPt",  &genPt,  "genPt/D");
    efficiencyTree->Branch("genEta", &genEta, "genEta/D");
    efficiencyTree->Branch("genPhi", &genPhi, "genPhi/D");

    // RCT cluster that was matched to the highest pT electron per event
    efficiencyTree->Branch("rct_cPt",  &rct_cPt,  "rct_cPt/D");
    efficiencyTree->Branch("rct_cEta", &rct_cEta, "rct_cEta/D");
    efficiencyTree->Branch("rct_cPhi", &rct_cPhi, "rct_cPhi/D");
    efficiencyTree->Branch("rct_deltaR", &rct_deltaR, "rct_deltaR/D");

    // GCT cluster that was matched to the highest pT electron per event
    efficiencyTree->Branch("gct_cPt",  &gct_cPt,  "gct_cPt/D");
    efficiencyTree->Branch("gct_cEta", &gct_cEta, "gct_cEta/D");
    efficiencyTree->Branch("gct_cPhi", &gct_cPhi, "gct_cPhi/D");
    efficiencyTree->Branch("gct_deltaR", &gct_deltaR, "gct_deltaR/D");
    
    nEvents = tfs_->make<TH1F>( "nEvents", "nEvents",  2,  0., 2. );

    // Histograms for rates
    l1eg_pt       = tfs_->make<TH1F>( "l1eg_pt"        , "p_{t}", 150,  0., 150. );
    l1egVLoose_pt = tfs_->make<TH1F>( "l1egVLoose_pt"  , "p_{t}", 150,  0., 150. );
    l1egLoose_pt  = tfs_->make<TH1F>( "l1egLoose_pt"   , "p_{t}", 150,  0., 150. );
    l1egMedium_pt = tfs_->make<TH1F>( "l1egMedium_pt"  , "p_{t}", 150,  0., 150. );
    l1egTight_pt  = tfs_->make<TH1F>( "l1egTight_pt"   , "p_{t}", 150,  0., 150. );

  }

void L1TCaloEGammaAnalyzerRates::beginJob( const EventSetup & es) {
}

void L1TCaloEGammaAnalyzerRates::analyze( const Event& evt, const EventSetup& es )
 {

  run = evt.id().run();
  lumi = evt.id().luminosityBlock();
  event = evt.id().event();
  
  nEvents->Fill(1);

  edm::Handle<l1tp2::CaloCrystalClusterCollection> rctCaloCrystalClusters;
  edm::Handle<l1tp2::CaloTowerCollection> rctCaloL1Towers;
  
  edm::Handle<l1tp2::CaloCrystalClusterCollection> gctCaloCrystalClusters;
  edm::Handle<l1tp2::CaloTowerCollection> gctCaloL1Towers;
  
  edm::Handle<EcalEBTrigPrimDigiCollection> ecalTPGs;
  edm::Handle<HcalTrigPrimDigiCollection> hcalTPGs;  
  edm::Handle<edm::SortedCollection<HcalTriggerPrimitiveDigi> > hbhecoll;
 
  std::vector<TLorentzVector> rctClustersMatched;
  std::vector<TLorentzVector> gctClustersMatched;

  rctClusters->clear(); 
  rctTowers->clear();
  gctClusters->clear();
  gctTowers->clear();
  allEcalTPGs->clear(); 
  allHcalTPGs->clear(); 

  // Detector geometry
  es.get<CaloGeometryRecord>().get(caloGeometry_);
  ebGeometry = caloGeometry_->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);
  hbGeometry = caloGeometry_->getSubdetectorGeometry(DetId::Hcal, HcalBarrel);
  es.get<HcalRecNumberingRecord>().get(hbTopology);
  hcTopology_ = hbTopology.product();
  HcalTrigTowerGeometry theTrigTowerGeometry(hcTopology_);
  es.get<CaloTPGRecord>().get(decoder_);

  std::cout << "Doing event " << event << "...." << std::endl;

  // Get the RCT clusters from the emulator and sort them by pT
  if(evt.getByToken(rctClustersSrc_, rctCaloCrystalClusters)){
    for(const auto & rctCluster : *rctCaloCrystalClusters){

      TLorentzVector temp ;
      std::cout << "RCT Cluster found: pT " << rctCluster.pt()  << ", "
		<< "eta "                   << rctCluster.eta() << ", "
		<< "phi "                   << rctCluster.phi() << std::endl;
      temp.SetPtEtaPhiE(rctCluster.pt(),rctCluster.eta(),rctCluster.phi(),rctCluster.pt());
      rctClusters->push_back(temp);
    }
  }
  std::sort(rctClusters->begin(), rctClusters->end(), L1TCaloEGammaAnalyzerRates::comparePt);

  // Get the RCT towers from the emulator 
  if(evt.getByToken(rctTowersSrc_, rctCaloL1Towers)) {
    for (const auto & rctTower : *rctCaloL1Towers){
      TLorentzVector temp;
      float totalEt = rctTower.ecalTowerEt() + rctTower.hcalTowerEt();
      if (totalEt > 0) {
	// std::cout << "Tower found: ECAL ET " << rctTower.ecalTowerEt()  << ", "
	// 	  << "HCAL ET " << rctTower.hcalTowerEt()  << ", "
	// 	  << "iEta, iPhi " << rctTower.towerIEta() << " " << rctTower.towerIPhi() << ", "
	// 	  << "eta "             << rctTower.towerEta() << ", "
	// 	  << "phi "             << rctTower.towerPhi() << std::endl;
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
      TLorentzVector temp ;
      std::cout << "GCT Cluster found: pT " << gctCluster.pt()  << ", "
                << "eta "               << gctCluster.eta() << ", "
                << "phi "               << gctCluster.phi() << std::endl;
      temp.SetPtEtaPhiE(gctCluster.pt(),gctCluster.eta(),gctCluster.phi(),gctCluster.pt());
      gctClusters->push_back(temp);
    }
  }
  std::sort(gctClusters->begin(), gctClusters->end(), L1TCaloEGammaAnalyzerRates::comparePt);

  // // get the ECAL inputs (i.e. ECAL crystals)
  // if(!evt.getByToken(ecalSrc_, ecalTPGs))
  //   std::cout<<"ERROR GETTING THE ECAL TPGS"<<std::endl;
  // else
  //   for (const auto& hit : *ecalTPGs.product()) {
  //     if (hit.encodedEt() > 0)  // hit.encodedEt() returns an int corresponding to 2x the crystal Et
  // 	{
  // 	  // Et is 10 bit, by keeping the ADC saturation Et at 120 GeV it means that you have to divide by 8
  // 	  float et = hit.encodedEt() / 8.;
	  
  // 	  if (et < 0.5)
  // 	    continue;  // keep the 500 MeV ET Cut
	  
  // 	  // std::cout << "ECAL hit et: " << et << std::endl;
	  
  // 	  // Get cell coordinates and info
  // 	  auto cell = ebGeometry->getGeometry(hit.id());
	  
  // 	  // std::cout << "Found ECAL cell/hit with coordinates " << cell->getPosition().x() << "," 
  // 	  // 	  << cell->getPosition().y() << "," 
  // 	  // 	  << cell->getPosition().z() << " and ET (GeV) " 
  // 	  // 	  << et << std::endl;
	  
  // 	  GlobalVector position=GlobalVector(cell->getPosition().x(), cell->getPosition().y(), cell->getPosition().z());
  // 	  float eta = position.eta();
  // 	  float phi = position.phi();
  // 	  TLorentzVector temp ;
  // 	  temp.SetPtEtaPhiE(et,eta,phi,et);
	  
  // 	  allEcalTPGs->push_back(temp);
  // 	}
  //   }
  
  // //ESHandle<L1CaloHcalScale> hcalScale;
  // //es.get<L1CaloHcalScaleRcd>().get(hcalScale);

  // if(!evt.getByToken(hcalSrc_, hcalTPGs))
  //   std::cout<<"ERROR GETTING THE HCAL TPGS"<<std::endl;
  // else
  // for (const auto& hit : *hcalTPGs.product()) {
  //   float et = decoder_->hcaletValue(hit.id(), hit.t0());
  //   ap_uint<10> encodedEt = hit.t0().compressedEt(); 
  //   // same thing as SOI_compressedEt() in HcalTriggerPrimitiveDigi.h///
  //   if (et <= 0)
  //     continue;
    
  //   if (!(hcTopology_->validHT(hit.id()))) {
  //     LogError("Phase2L1CaloEGammaEmulator")
  // 	<< " -- Hcal hit DetID not present in HCAL Geom: " << hit.id() << std::endl;
  //     throw cms::Exception("Phase2L1CaloEGammaEmulator");
  //     continue;
  //   }
  //   const std::vector<HcalDetId>& hcId = theTrigTowerGeometry.detIds(hit.id());
  //   if (hcId.empty()) {
  //     LogError("Phase2L1CaloEGammaEmulator")
  // 	<< "Cannot find any HCalDetId corresponding to " << hit.id() << std::endl;
  //     throw cms::Exception("Phase2L1CaloEGammaEmulator");
  //     continue;
  //   }
  //   if (hcId[0].subdetId() > 1)
  //     continue;
  //   GlobalVector hcal_tp_position = GlobalVector(0., 0., 0.);
  //   for (const auto& hcId_i : hcId) {
  //     if (hcId_i.subdetId() > 1)
  //       continue;
  //     // get the first HCAL TP/ cell
  //     auto cell = hbGeometry->getGeometry(hcId_i);
  //     if (cell == nullptr)
  // 	continue;
  //     GlobalVector tmpVector = GlobalVector(cell->getPosition().x(), cell->getPosition().y(), cell->getPosition().z());
  //     hcal_tp_position = tmpVector;
      
  //     // std::cout << "Found HCAL cell/TP with coordinates " << cell->getPosition().x() << ","
  //     //  		<< cell->getPosition().y() << ","
  //     //  		<< cell->getPosition().z() << " and ET (GeV) " << et
  //     // 		<< ", encoded Et " << encodedEt << std::endl;
      
  //     break;
  //   }
  
  //   float eta = hcal_tp_position.eta();
  //   float phi = hcal_tp_position.phi();
  //   TLorentzVector temp ;
  //   temp.SetPtEtaPhiE(et,eta,phi,et);
  //   allHcalTPGs->push_back(temp);
  // }


  //************************************************************************************/ 
  // Build the collections that we need for rates: only use the highest pT GCT cluster in each event
  //************************************************************************************/ 

  std::cout << "Building collections for rates: there are " << gctClusters->size() << " GCT clusters in the event" << std::endl;
  if (gctClusters->size() > 0) {
    
    l1eg_pt->Fill(gctClusters->at(0).Pt());
  
    if (gctClusters->at(0).Pt() > 25) {
      l1egVLoose_pt->Fill(gctClusters->at(0).Pt());
    }
    
    if (gctClusters->at(0).Pt() > 30) {
      l1egLoose_pt->Fill(gctClusters->at(0).Pt());
    }
    
    if (gctClusters->at(0).Pt() > 35) {
      l1egMedium_pt->Fill(gctClusters->at(0).Pt());
    }
    
    if (gctClusters->at(0).Pt() > 40) {
      l1egTight_pt->Fill(gctClusters->at(0).Pt());
    }
  }
  
 }



void L1TCaloEGammaAnalyzerRates::endJob() {
}

L1TCaloEGammaAnalyzerRates::~L1TCaloEGammaAnalyzerRates(){
}

DEFINE_FWK_MODULE(L1TCaloEGammaAnalyzerRates);

