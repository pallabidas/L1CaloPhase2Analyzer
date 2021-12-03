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


// ECAL propagation
#include "CommonTools/BaseParticlePropagator/interface/BaseParticlePropagator.h"
#include "CommonTools/BaseParticlePropagator/interface/RawParticle.h"


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

    // Gen electrons
    efficiencyTree->Branch("genPt",  &genPt,  "genPt/D");
    efficiencyTree->Branch("genEta", &genEta, "genEta/D");
    efficiencyTree->Branch("genPhi", &genPhi, "genPhi/D");

    // The RCT emulator cluster that was matched to the gen electron
    efficiencyTree->Branch("rct_cPt",  &rct_cPt,  "rct_cPt/D");
    efficiencyTree->Branch("rct_cEta", &rct_cEta, "rct_cEta/D");
    efficiencyTree->Branch("rct_cPhi", &rct_cPhi, "rct_cPhi/D");
    efficiencyTree->Branch("rct_deltaR", &rct_deltaR, "rct_deltaR/D");

    // The GCT emulator cluster that was matched to the gen electron
    efficiencyTree->Branch("gct_cPt",  &gct_cPt,  "gct_cPt/D");
    efficiencyTree->Branch("gct_cEta", &gct_cEta, "gct_cEta/D");
    efficiencyTree->Branch("gct_cPhi", &gct_cPhi, "gct_cPhi/D");
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
      // Not sure what this is in the code snippet
      // TLorentzVector trueElectron = TLorentzVector(prop.particle().E()*sin(prop.particle().vertex().theta()),
      // 						   prop.particle().vertex().eta(),
      // 						   prop.particle().vertex().phi(),
      // 						   0.);
      GlobalPoint ecal_pos(prop.particle().vertex().x(), prop.particle().vertex().y(), prop.particle().vertex().z());
      TLorentzVector corrGenElectron;
      // corrGenElectron.SetPtEtaPhiM(prop.particle().Pt(),
      // 				   prop.particle().eta(),
      // 				   prop.particle().phi(),
      // 				   prop.particle().mass());

      corrGenElectron.SetPtEtaPhiM(prop.particle().Pt(),
				   ecal_pos.eta(),
				   ecal_pos.phi(),
				   prop.particle().mass());
      
      std::cout << ">>> Corrected gen electron pt/eta/phi: (" 
		<< genElectron.pt() << ", " << genElectron.eta() << ", " << genElectron.phi() 
		<< ")"
		<< " to " 
		<< corrGenElectron.Pt() << ", " << corrGenElectron.Eta() << ", " << corrGenElectron.Phi()
		<< " with ecal_pos eta/phi "
		<< ecal_pos.eta() << ", " << ecal_pos.phi()
		<< std::endl;
      propagatedGenElectrons.push_back(corrGenElectron);
    }
  }

  //************************************************************************************/ 
  // Loop through the gen-level clusters and match to the RCT clusters, then the GCT clusters.
  //************************************************************************************/ 
  for (auto genElectron : propagatedGenElectrons) {

    genPt = genElectron.Pt();
    genEta = genElectron.Eta();
    genPhi = genElectron.Phi();

    rct_cPt    = 0;    rct_cEta   = -999;    rct_cPhi   = -999;
    rct_deltaR = 999;

    std::vector<TLorentzVector> rctClustersMatched;
    std::vector<TLorentzVector> gctClustersMatched;

    // Loop through the RCT clusters which are already sorted in decreasing pT, and check for
    // the first cluster within deltaR < 0.5. 
    std::cout << "Event " << event << ": check that RCT clusters are sorted!" << std::endl;
    for (size_t i = 0; i < rctClusters->size(); ++i) {
      std::cout << "sorted? RCT Cluster found: pT " << rctClusters->at(i).Pt()  << ", "
                << "eta "                   << rctClusters->at(i).Eta() << ", "
                << "phi "                   << rctClusters->at(i).Phi() << std::endl;
    }
    
    for (size_t i = 0; i < rctClusters->size(); ++i) {

      float this_rct_deltaR = reco::deltaR(rctClusters->at(i).Eta(), rctClusters->at(i).Phi(),
					   genElectron.Eta(), genElectron.Phi());
      // std::cout << "   Comparing "<< this_rct_deltaR << " to current rct_deltaR " << rct_deltaR << std::endl;

      if (this_rct_deltaR < 0.5) {

	// std::cout << "rctClusters pT " << rctClusters->at(i).Pt() 
	// 	  << " eta " << rctClusters->at(i).Eta()
	// 	  << " phi " << rctClusters->at(i).Phi() << std::endl;
	
	TLorentzVector temp ;
	temp.SetPtEtaPhiE(rctClusters->at(i).Pt(), rctClusters->at(i).Eta(),
			  rctClusters->at(i).Phi(), rctClusters->at(i).M());
	rctClustersMatched.push_back(temp);
      }
    }
    
    // For this gen electron, sort the matched clusters by pT, and only save the highest pT one
    std::sort(rctClustersMatched.begin(), rctClustersMatched.end(), L1TCaloEGammaAnalyzer::comparePt);
    rct_cPt  = rctClustersMatched.at(0).Pt();
    rct_cEta = rctClustersMatched.at(0).Eta();
    rct_cPhi = rctClustersMatched.at(0).Phi();
    rct_deltaR = reco::deltaR(rct_cEta, rct_cPhi,
			      genElectron.Eta(), genElectron.Phi());
    std::cout << "--> Matched RCT cluster " << rct_cPt 
	      << " eta: " << rct_cEta
	      << " phi: " << rct_cPhi
	      << " with genElectron " << genPt
	      << " eta: " << genEta
	      << " phi: " << genPhi << std::endl;
    
    //************************************************************************************/
    // Loop through the GCT cluster 4-vectors, and gen-match. Only take the first 
    // GCT cluster which has deltaR < 0.5.
    //************************************************************************************/ 
    std::cout << ">>> Doing GCT clusters..." << std::endl;
    
    gct_cPt    = 0;    gct_cEta   = -999;   gct_cPhi  = -999;
    gct_deltaR = 999;
    
    for (size_t i = 0; i < gctClusters->size(); ++i) {
      // std::cout << " gctClusters pT " << gctClusters->at(i).Pt() 
      // 	  << " eta " << gctClusters->at(i).Eta()
      // 	  << " phi " << gctClusters->at(i).Phi() << std::endl;
      float this_gct_deltaR = reco::deltaR(gctClusters->at(i).Eta(), gctClusters->at(i).Phi(),
					   genElectron.Eta(), genElectron.Phi());
      if (this_gct_deltaR < 0.5) {
	TLorentzVector temp;
        temp.SetPtEtaPhiE(gctClusters->at(i).Pt(), gctClusters->at(i).Eta(),
                          gctClusters->at(i).Phi(), gctClusters->at(i).M());
        gctClustersMatched.push_back(temp);
      }
    }
	
    // For this gen electron, sort the matched clusters by pT, and only save the highest pT one     
    std::sort(gctClustersMatched.begin(), gctClustersMatched.end(), L1TCaloEGammaAnalyzer::comparePt);
    gct_cPt  = gctClustersMatched.at(0).Pt();
    gct_cEta = gctClustersMatched.at(0).Eta();
    gct_cPhi = gctClustersMatched.at(0).Phi();
    gct_deltaR = reco::deltaR(gct_cEta, gct_cPhi,
                              genElectron.Eta(), genElectron.Phi());
    std::cout << "--> Matched GCT cluster " << gct_cPt
              << " eta: " << gct_cEta
              << " phi: " << gct_cPhi
              << " with genElectron " << genPt
              << " eta: " << genEta
              << " phi: " << genPhi << std::endl;

    efficiencyTree->Fill();

  } // end of loop over gen electrons

 }



void L1TCaloEGammaAnalyzer::endJob() {
}

L1TCaloEGammaAnalyzer::~L1TCaloEGammaAnalyzer(){
}

DEFINE_FWK_MODULE(L1TCaloEGammaAnalyzer);
