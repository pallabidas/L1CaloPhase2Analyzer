/*
 *  \file L1TEventDisplayGenerator.cc
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
#include "DataFormats/L1TCalorimeterPhase2/interface/CaloPFCluster.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"

#include "L1Trigger/L1CaloTrigger/interface/ParametricCalibration.h"
#include "L1Trigger/L1TCalorimeter/interface/CaloTools.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "L1Trigger/L1CaloPhase2Analyzer/interface/L1TEventDisplayGenerator.h"
#include "DataFormats/Math/interface/deltaR.h"

//#include "L1Trigger/L1CaloTrigger/plugins/Phase2L1CaloEGammaEmulator.h"

using namespace edm;
using std::cout;
using std::endl;
using std::vector;

L1TEventDisplayGenerator::L1TEventDisplayGenerator( const ParameterSet & cfg ) :
  decoderToken_(esConsumes<CaloTPGTranscoder, CaloTPGRecord>(edm::ESInputTag("", ""))),
  caloGeometryToken_(esConsumes<CaloGeometry, CaloGeometryRecord>(edm::ESInputTag("", ""))),
  hbTopologyToken_(esConsumes<HcalTopology, HcalRecNumberingRecord>(edm::ESInputTag("", ""))),
  ecalSrc_(consumes<EcalEBTrigPrimDigiCollection>(cfg.getParameter<edm::InputTag>("ecalDigis"))),
  hcalSrc_(consumes<HcalTrigPrimDigiCollection>(cfg.getParameter<edm::InputTag>("hcalDigis"))),
  ecalClustersSrc_(consumes<l1tp2::CaloCrystalClusterCollection >(cfg.getParameter<edm::InputTag>("clusters"))),
  caloTowersSrc_(consumes<l1tp2::CaloTowerCollection >(cfg.getParameter<edm::InputTag>("towers"))),
  caloPFClustersSrc_(consumes<l1tp2::CaloPFClusterCollection >(cfg.getParameter<edm::InputTag>("PFclusters"))),
  //genSrc_(consumes<reco::GenParticleCollection>(cfg.getParameter<edm::InputTag>( "GenParticles"))),  // does not work
  //genSrc_(consumes<pat::PackedGenParticleCollection>(cfg.getParameter<edm::InputTag>( "GenParticles"))),
  genSrc_(consumes<std::vector<reco::GenParticle>>(cfg.getParameter<edm::InputTag>( "GenParticles"))),
  tauSrc_(consumes<pat::TauCollection>(cfg.getParameter<edm::InputTag>( "Taus"))),
  jetSrc_(consumes<pat::JetCollection>(cfg.getParameter<edm::InputTag>( "Jets")))
  {
    folderName_          = cfg.getUntrackedParameter<std::string>("folderName");
    efficiencyTree = tfs_->make<TTree>("efficiencyTree", "Efficiency Tree");

    //efficiencyTree->Branch("sumTpgs_Pt",  &sumTpgs_Pt); 
    //efficiencyTree->Branch("sumTpgs_Eta", &sumTpgs_Eta); 
    //efficiencyTree->Branch("sumTpgs_Phi", &sumTpgs_Phi); 

    ////putting bufsize at 32000 and changing split level to 0 so that the branch isn't split into multiple branches
    efficiencyTree->Branch("ecalClusters", "vector<TLorentzVector>", &ecalClusters, 32000, 0); 
    efficiencyTree->Branch("caloTowers",   "vector<TLorentzVector>", &caloTowers, 32000, 0);
    efficiencyTree->Branch("caloPFClusters", "vector<TLorentzVector>", &caloPFClusters, 32000, 0);
    efficiencyTree->Branch("hcalTPGs", "vector<TLorentzVector>", &allHcalTPGs, 32000, 0); 
    efficiencyTree->Branch("ecalTPGs", "vector<TLorentzVector>", &allEcalTPGs, 32000, 0); 
    efficiencyTree->Branch("genParticles", "vector<TLorentzVector>", &genParticles, 32000, 0);
    efficiencyTree->Branch("recoTaus", "vector<TLorentzVector>", &recoTaus, 32000, 0);
    efficiencyTree->Branch("recoJets", "vector<TLorentzVector>", &recoJets, 32000, 0);
    
    //efficiencyTree->Branch("signalPFCands", "vector<TLorentzVector>", &signalPFCands, 32000, 0); 
    //efficiencyTree->Branch("l1Jets", "vector<TLorentzVector>", &l1Jets, 32000, 0); 
    //efficiencyTree->Branch("recoJets", "vector<TLorentzVector>", &recoJets, 32000, 0); 

    efficiencyTree->Branch("run",    &run,     "run/I");
    efficiencyTree->Branch("lumi",   &lumi,    "lumi/I");
    efficiencyTree->Branch("event",  &event,   "event/I");
    efficiencyTree->Branch("nvtx",   &nvtx,         "nvtx/I");

  }

void L1TEventDisplayGenerator::beginJob( const EventSetup & es) {
}

void L1TEventDisplayGenerator::analyze( const Event& evt, const EventSetup& es )
 {

  run = evt.id().run();
  lumi = evt.id().luminosityBlock();
  event = evt.id().event();

  edm::Handle<l1tp2::CaloCrystalClusterCollection> caloCrystalClusters;
  edm::Handle<l1tp2::CaloTowerCollection> caloL1Towers;
  edm::Handle<l1tp2::CaloPFClusterCollection> PFClusters;
  
  edm::Handle<EcalEBTrigPrimDigiCollection> ecalTPGs;
  edm::Handle<HcalTrigPrimDigiCollection> hcalTPGs;  
  edm::Handle<edm::SortedCollection<HcalTriggerPrimitiveDigi> > hbhecoll;
 
  ecalClusters->clear(); 
  caloTowers->clear();
  caloPFClusters->clear();
  allEcalTPGs->clear(); 
  allHcalTPGs->clear(); 
  genParticles->clear();
  recoTaus->clear();
  recoJets->clear();

  // Detector geometry
  caloGeometry_ = &es.getData(caloGeometryToken_);
  ebGeometry = caloGeometry_->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);
  hbGeometry = caloGeometry_->getSubdetectorGeometry(DetId::Hcal, HcalBarrel);
  hcTopology_ = &es.getData(hbTopologyToken_);
  HcalTrigTowerGeometry theTrigTowerGeometry(hcTopology_);
  decoder_ = &es.getData(decoderToken_);

  // Get the gen particles
  bool barrelParticles = false;
  edm::Handle<reco::GenParticleCollection> genParticleHandle;
  //edm::Handle<pat::PackedGenParticleCollection> genParticleHandle;
  if(!evt.getByToken(genSrc_, genParticleHandle)) std::cout<<"No gen Particles Found "<<std::endl;
  if(evt.getByToken(genSrc_, genParticleHandle)){
    for (unsigned int i = 0; i< genParticleHandle->size(); i++){
      edm::Ptr<reco::GenParticle> ptr(genParticleHandle, i);
      //edm::Ptr<pat::PackedGenParticle> ptr(genParticleHandle, i);
      std::cout<<"gen particle: "<<ptr->pdgId()<<"\t"<<ptr->pt()<<"\t"<<ptr->eta()<<"\t"<<ptr->phi()<<"\t"<<ptr->energy()<<std::endl;
      if(abs(ptr->eta()) < 1.4841) barrelParticles = true;
      TLorentzVector temp ;
      temp.SetPtEtaPhiE(ptr->pt(), ptr->eta(), ptr->phi(), ptr->energy());
      genParticles->push_back(temp);
    }
  }

  // Get the taus
  edm::Handle<pat::TauCollection> tauHandle;
  if(!evt.getByToken(tauSrc_, tauHandle)) std::cout<<"No Taus Found "<<std::endl;
  if(evt.getByToken(tauSrc_, tauHandle)){
    for (unsigned int i = 0; i< tauHandle->size(); i++){
      edm::Ptr<pat::Tau> ptr(tauHandle, i);
      std::cout<<"tau: "<<ptr->pt()<<"\t"<<ptr->eta()<<"\t"<<ptr->phi()<<"\t"<<ptr->energy()<<std::endl;
      TLorentzVector temp ;
      temp.SetPtEtaPhiE(ptr->pt(), ptr->eta(), ptr->phi(), ptr->energy());
      recoTaus->push_back(temp);
    }
  }

  // Get the jets
  edm::Handle<pat::JetCollection> jetHandle;
  if(!evt.getByToken(jetSrc_, jetHandle)) std::cout<<"No Jets Found "<<std::endl;
  if(evt.getByToken(jetSrc_, jetHandle)){
    for (unsigned int i = 0; i< jetHandle->size(); i++){
      edm::Ptr<pat::Jet> ptr(jetHandle, i);
      std::cout<<"jet: "<<ptr->pt()<<"\t"<<ptr->eta()<<"\t"<<ptr->phi()<<"\t"<<ptr->energy()<<std::endl;
      TLorentzVector temp ;
      temp.SetPtEtaPhiE(ptr->pt(), ptr->eta(), ptr->phi(), ptr->energy());
      recoJets->push_back(temp);
    }
  }

  // Get the ECAL clusters
  if(evt.getByToken(ecalClustersSrc_, caloCrystalClusters)){
    for(const auto & caloCluster : *caloCrystalClusters){
      //    for( vector<l1tp2::CaloCrystalClusterCollection>::const_iterator caloCluster = caloCrystalClusters->begin(); 
      //	 caloCluster != caloCrystalClusters->end(); 
      //	 caloCluster++ ) {
      //fill vector
      TLorentzVector temp ;
//      std::cout << "Cluster found: pT " << caloCluster.pt()  << ", "
//		<< "eta "               << caloCluster.eta() << ", "
//		<< "phi "               << caloCluster.phi() << std::endl;
      temp.SetPtEtaPhiE(caloCluster.pt(),caloCluster.eta(),caloCluster.phi(),caloCluster.pt());
      ecalClusters->push_back(temp);
    }
  }

  // Get the Calo towers
  if(evt.getByToken(caloTowersSrc_, caloL1Towers)) {
    for (const auto & caloTower : *caloL1Towers){
      TLorentzVector temp;
      float totalEt = caloTower.ecalTowerEt() + caloTower.hcalTowerEt();
      if (totalEt > 0) {
//	std::cout << "Tower found: ECAL ET " << caloTower.ecalTowerEt()  << ", "
//		  << "HCAL ET " << caloTower.hcalTowerEt()  << ", "
//		  << "iEta, iPhi " << caloTower.towerIEta() << " " << caloTower.towerIPhi() << ", "
//		  << "eta "             << caloTower.towerEta() << ", "
//		  << "phi "             << caloTower.towerPhi() << std::endl;
	temp.SetPtEtaPhiE(totalEt,
			  caloTower.towerEta(), caloTower.towerPhi(),
			  totalEt);
	caloTowers->push_back(temp);
      }
    }
  }

  // Get PF clusters
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
  /*
    for (size_t i = 0; i < ecalTPGs->size(); ++i) {
      int cal_ieta = (*ecalTPGs)[i].id().ieta();
      int cal_iphi = (*ecalTPGs)[i].id().iphi();
      if(cal_iphi==0)
	std::cout<<"cal_phi is 0"<<std::endl;
      if(cal_ieta<-28)
	continue;
      if(cal_ieta>28)
	continue;
      int ieta = TPGEtaRange(cal_ieta);
      short zside = (*ecalTPGs)[i].id().zside();
      // TPG iPhi starts at 1 and goes to 72.  Let's index starting at zero.
      // TPG ieta ideal goes from 0-55.
      double LSB = 0.5;
      //      double et= (*ecalTPGs)[i].compressedEt()*LSB;
      double et= (*ecalTPGs)[i].encodedEt()*LSB;
      if(ieta<0){
	std::cout<<"sorry, ieta less than 1 :("<<std::endl;
	std::cout<<"cal_ieta "<<cal_ieta<<" ieta "<<ieta<<std::endl;
      }
      float eta = getRecoEta(ieta, zside);
      float phi = getRecoPhi(cal_iphi);
      if(et==0)
	continue;
      //std::cout<<"et "<<et<<std::endl;
      TLorentzVector temp ;
      temp.SetPtEtaPhiE(et,eta,phi,et);
      //if(et>5)
      //std::cout<<"Event Display tpg ecal pt() "<<temp.Pt()<< " eta " <<eta << " phi "<< phi <<std::endl;
      allEcalTPGs->push_back(temp);
      }*/

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
    /*
    for (size_t i = 0; i < hcalTPGs->size(); ++i) {
      HcalTriggerPrimitiveDigi tpg = (*hcalTPGs)[i];
      int cal_ieta = tpg.id().ieta();
      int cal_iphi = tpg.id().iphi();
      if(cal_ieta>28)continue; 
      if(cal_ieta<-28)continue; 
      int ieta = TPGEtaRange(cal_ieta);
      short absieta = std::abs(tpg.id().ieta());
      short zside = tpg.id().zside();
      double et = hcalScale->et(tpg.SOI_compressedEt(), absieta, zside); 
      //if(et>0)
      //std::cout<<"HCAL ET "<<et<<std::endl;
      if(ieta<0){
	std::cout<<"sorry, ieta less than 1 :("<<std::endl;
	std::cout<<"cal_ieta "<<cal_ieta<<" ieta "<<ieta<<std::endl;
      }
      float eta = getRecoEta(ieta, zside);
      float phi = getRecoPhi(cal_iphi);    
      TLorentzVector temp ;
      temp.SetPtEtaPhiE(et,eta,phi,et);
      allHcalTPGs->push_back(temp);
      }
    */
    float eta = hcal_tp_position.eta();
    float phi = hcal_tp_position.phi();
    TLorentzVector temp ;
    temp.SetPtEtaPhiE(et,eta,phi,et);
    allHcalTPGs->push_back(temp);
  }
  if(barrelParticles) efficiencyTree->Fill();
 }

int L1TEventDisplayGenerator::get5x5TPGs(const int maxTPGPt_eta, 
				     const int maxTPGPt_phi, 
				     const double eTowerETMap[73][57], 
				     const double hTowerETMap[73][57], 
				     std::vector<double>* hcalTpgs_pt, 
				     std::vector<double>* hcalTpgs_eta, 
				     std::vector<double>* hcalTpgs_phi, 
				     std::vector<double>* ecalTpgs_pt, 
				     std::vector<double>* ecalTpgs_eta, 
				     std::vector<double>* ecalTpgs_phi,
				     std::vector<double>* sumTpgs_pt, 
				     std::vector<double>* sumTpgs_eta, 
				     std::vector<double>* sumTpgs_phi){
  for (int j = -5; j < 6; ++j) {//phi
    for (int k = -5; k < 6; ++k) { //eta
      int tpgsquarephi= maxTPGPt_phi+j;
      int tpgsquareeta= maxTPGPt_eta+k;
      if (tpgsquarephi==-1) {tpgsquarephi=71;}
      if (tpgsquarephi==-2) {tpgsquarephi=70;}
      if (tpgsquarephi==-3) {tpgsquarephi=69;}
      if (tpgsquarephi==-4) {tpgsquarephi=68;}
      if (tpgsquarephi==-5) {tpgsquarephi=67;}
      if (tpgsquarephi==72) {tpgsquarephi=0;}
      if (tpgsquarephi==73) {tpgsquarephi=1;}
      if (tpgsquarephi==74) {tpgsquarephi=2;}
      if (tpgsquarephi==75) {tpgsquarephi=3;}
      if (tpgsquarephi==76) {tpgsquarephi=4;}
      if (tpgsquareeta>55 || tpgsquareeta<0) {continue;}//No Eta values beyond
      hcalTpgs_pt->push_back(hTowerETMap[tpgsquarephi][tpgsquareeta]);
      hcalTpgs_eta->push_back(towerEtaMap[k]);
      hcalTpgs_phi->push_back(towerPhiMap[j]);

      ecalTpgs_pt->push_back(eTowerETMap[tpgsquarephi][tpgsquareeta]);
      ecalTpgs_eta->push_back(towerEtaMap[tpgsquareeta]);
      ecalTpgs_phi->push_back(towerPhiMap[tpgsquarephi]);

      sumTpgs_pt->push_back(eTowerETMap[tpgsquarephi][tpgsquareeta]+eTowerETMap[tpgsquarephi][tpgsquareeta]);
      sumTpgs_eta->push_back(towerEtaMap[tpgsquareeta]);
      sumTpgs_phi->push_back(towerPhiMap[tpgsquarephi]);
    }
  }
  //std::cout<<"TPGe5x5_ "<<TPGe5x5_<<" TPGh5x5_ "<<TPGh5x5_<<std::endl;
  //return (TPGe5x5_ + TPGh5x5_);
  return 1;
}

/*
 * Get the ECAL TPGS create a TPG map for the event
 *
 */

void L1TEventDisplayGenerator::initializeECALTPGMap(Handle<EcalTrigPrimDigiCollection> ecal, double eTowerETMap[73][57], bool testMode){
  //std::cout << "ECAL TPGS" << std::endl;
  for (size_t i = 0; i < ecal->size(); ++i) {
    int cal_ieta = (*ecal)[i].id().ieta();
    int cal_iphi = (*ecal)[i].id().iphi();
    int iphi = cal_iphi-1;
    int ieta = TPGEtaRange(cal_ieta);
    // TPG iPhi starts at 1 and goes to 72.  Let's index starting at zero.
    // TPG ieta ideal goes from 0-55.
    double LSB = 0.5;
    double et= (*ecal)[i].compressedEt()*LSB;

    if(testMode && iphi == 34 && ieta == 11){
      et = 40;
    }

    if (iphi >= 0 && iphi <= 72 &&
	ieta >= 0 && ieta <= 55) {
      eTowerETMap[iphi][ieta] = et; 
    }

  }

}

void L1TEventDisplayGenerator::initializeHCALTPGMap(const Handle<HcalTrigPrimDigiCollection> hcal, 
					 const ESHandle<L1CaloHcalScale>  hcalScale, 
					 double hTowerETMap[73][57], bool testMode){
  for (size_t i = 0; i < hcal->size(); ++i) {
    HcalTriggerPrimitiveDigi tpg = (*hcal)[i];
    int cal_ieta = tpg.id().ieta();
    int cal_iphi = tpg.id().iphi();
    int iphi = cal_iphi-1;
    int ieta = TPGEtaRange(cal_ieta);
    short absieta = std::abs(tpg.id().ieta());
    short zside = tpg.id().zside();
    double energy = hcalScale->et(tpg.SOI_compressedEt(), absieta, zside); 

    if(testMode && iphi == 34 && ieta == 12){
      energy = 40;
    }

    if (iphi >= 0 && iphi <= 71 &&
	ieta >= 0 && ieta <= 55) {
      //(*hcal)[i].SOI_compressedEt(), absieta, zside)*LSB; //*LSB
      //if(energy>0)
      //std::cout<<"hcal iphi "<<iphi<<" ieta "<<ieta<<" energy "<<energy<<std::endl;
      hTowerETMap[iphi][ieta] = energy;
      //TPGSum_ +=energy;
      //TPGH_ += energy;
      //double alpha_h = TPGSFp_[cal_ieta]; //v3
      //hCorrTowerETMap[cal_iphi][cal_ieta] = alpha_h*energy;
      //cTPGH_ += alpha_h*energy;
      //if (energy > 0) {
      //std::cout << "hcal eta/phi=" << ieta << "/" << iphi
      //<< " = (" << getEtaTPG(ieta) << "/" << getPhiTPG(iphi) << ") "
      //<< " et=" << (*hcal)[i].SOI_compressedEt()
      //<< " energy=" << energy
      //<< " rctEta="<< twrEta2RegionEta(cal_ieta) << " rctPhi=" << twrPhi2RegionPhi(cal_iphi)
      //<< " fg=" << (*hcal)[i].SOI_fineGrain() << std::endl;
      //}
      //if (energy>maxTPGHPt){
      //maxTPGHPt=energy;
      //maxTPGHPt_phi = cal_iphi; //this one starts at 0-72
      //maxTPGHPt_eta = cal_ieta; //this one is 0-54
      //}
    }
    //else
      //std::cout<<"HCAL failed checks iphi "<<iphi<<" ieta "<<ieta<<std::endl;
  }//end HCAL TPG
}

void L1TEventDisplayGenerator::endJob() {
}

L1TEventDisplayGenerator::~L1TEventDisplayGenerator(){
}

DEFINE_FWK_MODULE(L1TEventDisplayGenerator);
