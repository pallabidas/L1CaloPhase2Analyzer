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

#include "L1Trigger/L1CaloPhase2Analyzer/interface/L1TCaloEGammaAnalyzer.h"
#include "DataFormats/Math/interface/deltaR.h"

//#include "L1Trigger/L1CaloTrigger/plugins/Phase2L1CaloEGammaEmulator.h"

using namespace edm;
using std::cout;
using std::endl;
using std::vector;

L1TCaloEGammaAnalyzer::L1TCaloEGammaAnalyzer( const ParameterSet & cfg ) :
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
    efficiencyTree->Branch("nvtx",   &nvtx,         "nvtx/I");

    // RCT emulator clusters and their matched gen electrons
    efficiencyTree->Branch("rct_cPt",  &rct_cPt,  "rct_cPt/D");
    efficiencyTree->Branch("rct_cEta", &rct_cEta, "rct_cEta/D");
    efficiencyTree->Branch("rct_cPhi", &rct_cPhi, "rct_cPhi/D");
    
    efficiencyTree->Branch("rct_genPt",  &rct_genPt,  "rct_genPt/D");
    efficiencyTree->Branch("rct_genEta", &rct_genEta, "rct_genEta/D");
    efficiencyTree->Branch("rct_genPhi", &rct_genPhi, "rct_genPhi/D");

    efficiencyTree->Branch("rct_deltaR", &rct_deltaR, "rct_deltaR/D");

    // GCT emulator clusters etc.
    efficiencyTree->Branch("gct_cPt",  &gct_cPt,  "gct_cPt/D");
    efficiencyTree->Branch("gct_cEta", &gct_cEta, "gct_cEta/D");
    efficiencyTree->Branch("gct_cPhi", &gct_cPhi, "gct_cPhi/D");

    efficiencyTree->Branch("gct_genPt",  &gct_genPt,  "gct_genPt/D");
    efficiencyTree->Branch("gct_genEta", &gct_genEta, "gct_genEta/D");
    efficiencyTree->Branch("gct_genPhi", &gct_genPhi, "gct_genPhi/D");

    efficiencyTree->Branch("gct_deltaR", &gct_deltaR, "gct_deltaR/D");
    
  }

void L1TCaloEGammaAnalyzer::beginJob( const EventSetup & es) {
}

void L1TCaloEGammaAnalyzer::analyze( const Event& evt, const EventSetup& es )
 {

  run = evt.id().run();
  lumi = evt.id().luminosityBlock();
  event = evt.id().event();

  edm::Handle<l1tp2::CaloCrystalClusterCollection> rctCaloCrystalClusters;
  edm::Handle<l1tp2::CaloTowerCollection> rctCaloL1Towers;
  
  edm::Handle<l1tp2::CaloCrystalClusterCollection> gctCaloCrystalClusters;
  edm::Handle<l1tp2::CaloTowerCollection> gctCaloL1Towers;
  
  edm::Handle<EcalEBTrigPrimDigiCollection> ecalTPGs;
  edm::Handle<HcalTrigPrimDigiCollection> hcalTPGs;  
  edm::Handle<edm::SortedCollection<HcalTriggerPrimitiveDigi> > hbhecoll;
 
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

  // Get the RCT clusters from the emulator and sort them by pT
  if(evt.getByToken(rctClustersSrc_, rctCaloCrystalClusters)){
    for(const auto & rctCluster : *rctCaloCrystalClusters){
      //    for( vector<l1tp2::CaloCrystalClusterCollection>::const_iterator rctCluster = caloCrystalClusters->begin(); 
      //	 rctCluster != caloCrystalClusters->end(); 
      //	 rctCluster++ ) {
      //fill vector
      TLorentzVector temp ;
      std::cout << "RCT Cluster found: pT " << rctCluster.pt()  << ", "
		<< "eta "                   << rctCluster.eta() << ", "
		<< "phi "                   << rctCluster.phi() << std::endl;
      temp.SetPtEtaPhiE(rctCluster.pt(),rctCluster.eta(),rctCluster.phi(),rctCluster.pt());
      rctClusters->push_back(temp);
    }
  }
  std::sort(rctClusters->begin(), rctClusters->end(), L1TCaloEGammaAnalyzer::comparePt);

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
  std::sort(gctClusters->begin(), gctClusters->end(), L1TCaloEGammaAnalyzer::comparePt);

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
	  
	  // std::cout << "ECAL hit et: " << et << std::endl;
	  
	  // Get cell coordinates and info
	  auto cell = ebGeometry->getGeometry(hit.id());
	  
	  // std::cout << "Found ECAL cell/hit with coordinates " << cell->getPosition().x() << "," 
	  // 	  << cell->getPosition().y() << "," 
	  // 	  << cell->getPosition().z() << " and ET (GeV) " 
	  // 	  << et << std::endl;
	  
	  GlobalVector position=GlobalVector(cell->getPosition().x(), cell->getPosition().y(), cell->getPosition().z());
	  float eta = position.eta();
	  float phi = position.phi();
	  TLorentzVector temp ;
	  temp.SetPtEtaPhiE(et,eta,phi,et);
	  
	  // int cc = 28;
	  // if (getCrystal_phiID(position.phi()) <= getPhiMax_card(cc) &&
	  //     getCrystal_phiID(position.phi()) >= getPhiMin_card(cc) &&
	  //     getCrystal_etaID(position.eta()) <= getEtaMax_card(cc) &&
	  //     getCrystal_etaID(position.eta()) >= getEtaMin_card(cc)){

	  //   // TEMP: only use ECAL hits in Card 28                                                                             
	  //   allEcalTPGs->push_back(temp); 

	  //   if ((getCrystal_etaID(position.eta()) > 29) && (getCrystal_etaID(position.eta()) < 35)) {
	  //     std::cout << "[CARD " << cc << "]: Found ECAL cell/hit with eta/phi "
	  // 		<< position.eta() << ", "
	  // 		<< position.phi() << ", and in-detector phiID and etaID "
	  // 		<< getCrystal_phiID(position.phi()) << ", "
	  // 		<< getCrystal_etaID(position.eta()) << ", and ET (GeV) "
	  // 		<< et << std::endl;
	  //   }
	  // }
	  // TEMP: Debugging purposes only: only add ECAL hits in Card 28
	  allEcalTPGs->push_back(temp);
	}
    }
  
  //ESHandle<L1CaloHcalScale> hcalScale;
  //es.get<L1CaloHcalScaleRcd>().get(hcalScale);

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
      
      // std::cout << "Found HCAL cell/TP with coordinates " << cell->getPosition().x() << ","
      //  		<< cell->getPosition().y() << ","
      //  		<< cell->getPosition().z() << " and ET (GeV) " << et
      // 		<< ", encoded Et " << encodedEt << std::endl;
      
      break;
    }
  
    float eta = hcal_tp_position.eta();
    float phi = hcal_tp_position.phi();
    TLorentzVector temp ;
    temp.SetPtEtaPhiE(et,eta,phi,et);
    allHcalTPGs->push_back(temp);
  }

  // Get genParticles and build a vector of genElectrons
  edm::Handle<GenParticleCollectionType> genParticleHandle;
  if(!evt.getByToken(genToken_,genParticleHandle))
    std::cout<<"No gen Particles Found "<<std::endl;
  else
    std::cout<<"Gen Particles size "<<genParticleHandle->size()<<std::endl;
  
  std::vector<reco::GenParticle> genElectrons;
  std::vector<reco::GenParticle> genParticles;
  
  for (unsigned int i = 0; i< genParticleHandle->size(); i++){
    edm::Ptr<reco::GenParticle> ptr(genParticleHandle, i);
    genParticles.push_back(*ptr);
    
    // Get gen electrons in barrel + overlap
    if ( (abs(ptr->pdgId()) == 11) && ( abs(ptr->eta()) < 1.4841 )) {
      genElectrons.push_back(*ptr);
      std::cout << "Added genElectron " << ptr->pt() << std::endl;
    }
  }

  // Keep track of which genElectrons are already matched to a RCT cluster, and to a GCT cluster
  // 0 if false, 1 if true
  int isGenEleMatchedToRCT[genElectrons.size()] = {0};
  int isGenEleMatchedToGCT[genElectrons.size()] = {0};

  //************************************************************************************/ 
  // Loop through the RCT cluster 4-vectors, and gen-match. If multiple RCT clusters
  // satisfy the deltaR criterion for a single gen electron, only match the highest pT
  // cluster, and the rejected clusters will have genEta = -999 (and all default values for gen*).
  //************************************************************************************/ 
  for (size_t i = 0; i < rctClusters->size(); ++i) {

    rct_cPt    = 0;    rct_cEta   = -999;    rct_cPhi   = -999;
    rct_genPt  = 0;    rct_genEta = -999;    rct_genPhi = -999;
    rct_deltaR = 999;
    
    std::cout << "rctClusters pT " << rctClusters->at(i).Pt() 
	      << " eta " << rctClusters->at(i).Eta()
	      << " phi " << rctClusters->at(i).Phi() << std::endl;

    rct_cPt  = rctClusters->at(i).Pt();
    rct_cEta = rctClusters->at(i).Eta();
    rct_cPhi = rctClusters->at(i).Phi();

    int thisGenEle = 0;  // counter for current gen electron
    int lastGenEle = 0;  // the last gen electron to have been matched
    for (auto genElectron : genElectrons) {
      float this_rct_deltaR = reco::deltaR(rctClusters->at(i).Eta(), rctClusters->at(i).Phi(),
				      genElectron.eta(), genElectron.phi());
      std::cout << "   Comparing "<< this_rct_deltaR << " to current rct_deltaR " << rct_deltaR << std::endl;
      // If this rct_deltaR is the smallest one so far and this gen electron has not been already matched
      if ((this_rct_deltaR < rct_deltaR) && (isGenEleMatchedToRCT[thisGenEle] == 0)) {
	rct_genPt = genElectron.pt();
	rct_genEta = genElectron.eta();
	rct_genPhi = genElectron.phi();
	rct_deltaR = this_rct_deltaR;
	lastGenEle = thisGenEle;
      }
      // Increment our counter as we loop through the gen electrons
      thisGenEle += 1;
    }

    // Modify our record of whether this gen electron was matched (the record persists as we
    // loop through the RCT clusters)
    std::cout << "   Mark this gen electron as matched, set isGenEleMatchedToRCT to 1 at element "\
	      << lastGenEle << std::endl;
    isGenEleMatchedToRCT[lastGenEle] = 1;
      

    // n.b. if a cluster was never matched, its rct_cPt will be zero, rct_cEta will be 999 (all default)
    if (rct_cPt > 0) {
      std::cout << "--> Matched cluster " << rct_cPt 
		<< " eta: " << rct_cEta
		<< " phi: " << rct_cPhi
		<< " with genElectron " << rct_genPt
		<< " eta: " << rct_genEta
		<< " phi: " << rct_genPhi << std::endl;
    }
    else {
      std::cout << "--> No genElectron matched" << std::endl;
    }
      
    efficiencyTree->Fill();
  } // end of loop over RCT clusters

  
  //************************************************************************************/
  // Loop through the GCT cluster 4-vectors, and gen-match. If multiple GCT clusters
  // satisfy the deltaR criterion for a single gen electron, only match the highest pT
  // cluster, and the rejected clusters will have genEta = -999 (and all default values for gen*).
  //************************************************************************************/ 
  std::cout << ">>> Doing GCT clusters..." << std::endl;
  for (size_t i = 0; i < gctClusters->size(); ++i) {

    gct_cPt    = 0;    gct_cEta   = -999;   gct_cPhi  = -999;
    gct_genPt  = 0;    gct_genEta = -999;   gct_genPhi = -999;
    gct_deltaR = 999;
    
    std::cout << " gctClusters pT " << gctClusters->at(i).Pt() 
	      << " eta " << gctClusters->at(i).Eta()
	      << " phi " << gctClusters->at(i).Phi() << std::endl;

    gct_cPt  = gctClusters->at(i).Pt();
    gct_cEta = gctClusters->at(i).Eta();
    gct_cPhi = gctClusters->at(i).Phi();

    int thisGenEle = 0;  // counter for current gen electron
    int lastGenEle = 0;  // the last gen electron to have been matched
    for (auto genElectron : genElectrons) {
      float this_gct_deltaR = reco::deltaR(gctClusters->at(i).Eta(), gctClusters->at(i).Phi(),
					   genElectron.eta(), genElectron.phi());
      std::cout << "     Comparing "<< this_gct_deltaR << " to current gct_deltaR " << gct_deltaR << std::endl;
      // If this gct_deltaR is the smallest one so far and this gen electron has not been already matched
      if ((this_gct_deltaR < gct_deltaR) && (isGenEleMatchedToGCT[thisGenEle] == 0)) {
	gct_genPt = genElectron.pt();
	gct_genEta = genElectron.eta();
	gct_genPhi = genElectron.phi();
	gct_deltaR = this_gct_deltaR;
	lastGenEle = thisGenEle;
      }
      // Increment our counter as we loop through the gen electrons
      thisGenEle += 1;
    }

    // Modify our record of whether this gen electron was matched (the record persists as we
    // loop through the GCT clusters)
    std::cout << "   Mark this gen electron as matched, set isGenEleMatchedToGCT to 1 at element "\
	      << lastGenEle << std::endl;
    isGenEleMatchedToGCT[lastGenEle] = 1;
      

    // n.b. if a cluster was never matched, its gct_cPt will be zero, gct_cEta will be 999 (all default)
    if (gct_cPt > 0.0) {
      std::cout << "  --> Matched cluster " << gct_cPt 
		<< " eta: " << gct_cEta
		<< " phi: " << gct_cPhi
		<< " with genElectron " << gct_genPt
		<< " eta: " << gct_genEta
		<< " phi: " << gct_genPhi << std::endl;
    }
    else {
      std::cout << "--> No genElectron matched" << std::endl;
    }
      
    efficiencyTree->Fill();
  } // end of loop over GCT clusters
  
  

  
 }


void L1TCaloEGammaAnalyzer::endJob() {
}

L1TCaloEGammaAnalyzer::~L1TCaloEGammaAnalyzer(){
}

DEFINE_FWK_MODULE(L1TCaloEGammaAnalyzer);
