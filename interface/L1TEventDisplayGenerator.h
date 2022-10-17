#ifndef L1TEventDisplayGenerator_H
#define L1TEventDisplayGenerator_H


// system include files
#include <memory>
#include <unistd.h>


#include <iostream>
#include <fstream>
#include <vector>

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
//#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

// GCT and RCT data formats
#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"
#include "DataFormats/L1GlobalCaloTrigger/interface/L1GctCollections.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "DataFormats/TauReco/interface/BaseTau.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/PFTauTagInfo.h"

#include <memory>
#include <math.h>
#include <vector>
#include <list>
#include <TLorentzVector.h>

#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "CondFormats/L1TObjects/interface/L1CaloHcalScale.h"
#include "CondFormats/DataRecord/interface/L1CaloHcalScaleRcd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "L1Trigger/L1TGlobal/interface/TriggerMenuFwd.h"
#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloRegion.h"
#include "L1Trigger/L1TCaloLayer1/src/UCTRegion.hh"
#include "L1Trigger/L1TCaloLayer1/src/UCTGeometry.hh"
#include "DataFormats/L1TCalorimeterPhase2/interface/CaloCrystalCluster.h"

//
// class declaration
//
using std::vector;

static constexpr int n_crystals_towerEta = 5;
static constexpr int n_crystals_towerPhi = 5;
static constexpr int n_towers_Eta = 34;
static constexpr int n_towers_Phi = 72;
static constexpr float ECAL_eta_range = 1.4841;
static constexpr float half_crystal_size = 0.00873;
static constexpr int n_towers_cardEta = 17;   // new: equivalent to n_towers_per_link
static constexpr int n_towers_cardPhi = 4;    
static constexpr int n_towers_per_link = 17;

static constexpr int CRYSTALS_IN_TOWER_ETA = 5;
static constexpr int CRYSTALS_IN_TOWER_PHI = 5;

static constexpr int TOWER_IN_ETA = 3;      // number of towers in eta, in one 3x4 region (barrel)
static constexpr int TOWER_IN_PHI = 4;      // number of towers in phi, in one 3x4 region (barrel)

//static constexpr int TOWER_IN_ETA_OVERLAP = 2; // number of towers in eta, in one 2x4 region (overlap)
//static constexpr int TOWER_IN_PHI_OVERLAP = 4; // number of towers in phi, in one 2x4 region (overlap)

static constexpr int CRYSTAL_IN_ETA = 15;   // number of crystals in eta, in one 3x4 region (barrel)
static constexpr int CRYSTAL_IN_PHI = 20;   // number of crystals in phi, in one 3x4 region (barrel)


static constexpr int N_CLUSTERS_PER_REGION = 4;       // number of clusters per ECAL region
static constexpr int N_REGIONS_PER_CARD = 6;          // number of ECAL regions per card

// Assert that the card index is within bounds. (Valid cc: 0 to 35, since there are 36 RCT cards)
bool isValidCard(int cc) {
  return ((cc > -1) && (cc < 36));
}


int getCrystal_etaID(float eta) {
  float size_cell = 2 * ECAL_eta_range / (n_crystals_towerEta * n_towers_Eta);
  int etaID = int((eta + ECAL_eta_range) / size_cell);
  return etaID;
}

int getCrystal_phiID(float phi) {
  float size_cell = 2 * M_PI / (n_crystals_towerPhi * n_towers_Phi);
  int phiID = int((phi + M_PI) / size_cell);
  return phiID;
}


int getEtaMax_card(int card) {
  int etamax = 0;
  if (card % 2 == 0)
    etamax = n_towers_per_link * n_crystals_towerEta - 1;  // First eta half. 5 crystals in eta in 1 tower.
  else
    etamax = n_towers_Eta * n_crystals_towerEta - 1;
  return etamax;
}

int getEtaMin_card(int card) {
  int etamin = 0;
  if (card % 2 == 0)
    etamin = 0 * n_crystals_towerEta;  // First eta half. 5 crystals in eta in 1 tower.
  else
    etamin = n_towers_per_link * n_crystals_towerEta;
  return etamin;
}

int getPhiMax_card(int card) {
  int phimax = ((card / 2) + 1) * 4 * n_crystals_towerPhi - 1;
  return phimax;
}

int getPhiMin_card(int card) {
  int phimin = (card / 2) * 4 * n_crystals_towerPhi;
  return phimin;
}

class L1TEventDisplayGenerator : public edm::one::EDAnalyzer<edm::one::SharedResources> {

 public:
  
  // Constructor
  L1TEventDisplayGenerator(const edm::ParameterSet& ps);
  
  // Destructor
  virtual ~L1TEventDisplayGenerator();

  edm::Service<TFileService> tfs_;

  std::vector<double> *hcalTpgs_Pt  = new std::vector<double>; 
  std::vector<double> *hcalTpgs_Eta = new std::vector<double>; 
  std::vector<double> *hcalTpgs_Phi = new std::vector<double>; 

  std::vector<double> *ecalTpgs_Pt  = new std::vector<double>; 
  std::vector<double> *ecalTpgs_Eta = new std::vector<double>; 
  std::vector<double> *ecalTpgs_Phi = new std::vector<double>; 

  std::vector<double> *sumTpgs_Pt  = new std::vector<double>; 
  std::vector<double> *sumTpgs_Eta = new std::vector<double>; 
  std::vector<double> *sumTpgs_Phi = new std::vector<double>; 

  std::vector<TLorentzVector> *rlxTaus  = new std::vector<TLorentzVector>; 
  std::vector<TLorentzVector> *isoTaus  = new std::vector<TLorentzVector>; 
  std::vector<TLorentzVector> *recoTaus  = new std::vector<TLorentzVector>; 
  std::vector<TLorentzVector> *allRegions  = new std::vector<TLorentzVector>; 
  std::vector<TLorentzVector> *allEcalTPGs  = new std::vector<TLorentzVector>; 
  std::vector<TLorentzVector> *allHcalTPGs  = new std::vector<TLorentzVector>; 
  std::vector<TLorentzVector> *signalPFCands  = new std::vector<TLorentzVector>; 
  std::vector<TLorentzVector> *l1Jets  = new std::vector<TLorentzVector>; 
  std::vector<TLorentzVector> *recoJets  = new std::vector<TLorentzVector>; 
  std::vector<double> *recoJetsDr  = new std::vector<double>;

  // Outputs of the emulator
  std::vector<TLorentzVector> *ecalClusters  = new std::vector<TLorentzVector>; 
  std::vector<TLorentzVector> *caloTowers    = new std::vector<TLorentzVector>;
  std::vector<TLorentzVector> *caloPFClusters = new std::vector<TLorentzVector>;

  TH1F* isoTau_pt;
  TH1F* isoTau_eta;
  TH1F* isoTau_phi;

  TH1F* tau_pt;
  TH1F* tau_eta;
  TH1F* tau_phi;

  TH1F* recoTau_pt;
  TH1F* recoTau_eta;
  TH1F* recoTau_phi;
  TTree* efficiencyTree;

  int run, lumi, event;
  double isoTauPt, rlxTauPt, isoTauEta, rlxTauEta, isoTauPhi, rlxTauPhi;
  double recoPt, recoEta, recoPhi;
  int l1RlxMatched, l1IsoMatched;
  int decayMode;
  double tauEtaEcalEnt,tauPhiEcalEnt,rawEcal, rawHcal, ecal, hcal, jetEt, jetEta, jetPhi, nvtx;
  double max3ProngDeltaR, minProngPt, maxProngPt, midProngPt; int n3ProngCands;
  double pfCandsEt, signalCandsEt, isoCandsEt;
  double TPG2x2, TPGH2x2, TPGE2x2;
  double TPG5x5, TPGH5x5, TPGE5x5;
  double TPG6x6, TPGH6x6, TPGE6x6;
  double TPG7x7, TPGH7x7, TPGE7x7;

  void getThreeProngInfo(const pat::Tau & tau, double &maxDeltaR, double &minProngPt, double &midProngPt, double &maxProngPt, int &nCands);
  void getRawEcalHcalEnergy(const pat::PackedCandidate pfCand, double &rawEcal, double &rawHcal, double &ecal, double &hcal);
  double getPFCandsEt(const std::vector<pat::PackedCandidate> pfCands);
  double getPFCandsEtEtaPhi(edm::Handle<std::vector<pat::PackedCandidate> >& pfCands, const pat::Tau &tau, double dR);
  void initializeHCALTPGMap(const edm::Handle<HcalTrigPrimDigiCollection> hcal, const  edm::ESHandle<L1CaloHcalScale> hcalScale, double hTowerETMap[73][57], bool testMode = false);
  void initializeECALTPGMap(edm::Handle<EcalTrigPrimDigiCollection> ecal, double eTowerETMap[73][57], bool testMode = false);

int get5x5TPGs(const int maxTPGPt_eta, 
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
	       std::vector<double>* sumTpgs_phi);

 protected:
  // Analyze
  void analyze(const edm::Event& evt, const edm::EventSetup& es);
  
  // BeginJob
  void beginJob(const edm::EventSetup &es);
  
  // EndJob
  void endJob(void);

  
 private:
  // ----------member data ---------------------------

  int nev_; // Number of events processed
  bool verbose_;
  std::ofstream logFile_;
  edm::InputTag rctSource_; 

  edm::ESGetToken<CaloTPGTranscoder, CaloTPGRecord> decoderToken_;
  const CaloTPGTranscoder* decoder_;

  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeometryToken_;
  const CaloGeometry* caloGeometry_;
  const CaloSubdetectorGeometry* ebGeometry;
  const CaloSubdetectorGeometry* hbGeometry;
  edm::ESGetToken<HcalTopology, HcalRecNumberingRecord> hbTopologyToken_;
  const HcalTopology* hcTopology_;

  edm::EDGetTokenT<vector<pat::PackedCandidate> > packedPfCandsToken_;  
  edm::EDGetTokenT<vector<reco::PFCandidate> > pfCandsToken_;  
  edm::EDGetTokenT<L1CaloRegionCollection> L1RegionCollection;
  edm::EDGetTokenT<L1CaloEmCollection> L1EMCollection_;
  edm::EDGetTokenT<reco::VertexCollection> vertices_;
  edm::EDGetTokenT<EcalEBTrigPrimDigiCollection> ecalSrc_; 
  edm::EDGetTokenT<HcalTrigPrimDigiCollection> hcalSrc_;
  //edm::EDGetTokenT<double> recoPt_;
  //edm::EDGetTokenT<std::string> folderName_;
  edm::EDGetTokenT<reco::VertexCollection> vtxLabel_;
  edm::EDGetTokenT<reco::PFTauDiscriminator> discriminatorMu_;
  edm::EDGetTokenT<reco::PFTauDiscriminator> discriminatorIso_;
  edm::EDGetTokenT<vector<reco::PFTau> > tauSrc_;
  edm::EDGetTokenT<vector<pat::Tau> > slimmedTauSrc_;
  edm::EDGetTokenT<vector<pat::Jet> > jetSrc_;
  edm::EDGetTokenT<L1GctJetCandCollection> gctIsoTauJetsSource_;
  edm::EDGetTokenT<L1GctJetCandCollection> gctTauJetsSource_;
  edm::EDGetTokenT<vector <l1extra::L1JetParticle> > l1ExtraIsoTauSource_;
  edm::EDGetTokenT<vector <l1extra::L1JetParticle> > l1ExtraTauSource_;
  edm::EDGetTokenT<vector <l1extra::L1JetParticle> > l1ExtraJetSource_;
  std::vector< edm::EDGetTokenT<l1t::TauBxCollection> > stage2TauSource_;
  edm::EDGetTokenT<vector <L1CaloRegion> > regionSource_;
  edm::EDGetTokenT<l1tp2::CaloCrystalClusterCollection> ecalClustersSrc_;
  edm::EDGetTokenT<l1tp2::CaloTowerCollection> caloTowersSrc_;
  edm::EDGetTokenT<l1tp2::CaloPFClusterCollection> caloPFClustersSrc_;

  std::string folderName_;
  double recoPt_;

		 
  int TPGEtaRange(int ieta){
    int iEta = 0;
    // So here, -28 becomes 0.  -1 be comes 27.  +1 becomes 28. +28 becomes 55.
    // And we have mapped [-28, -1], [1, 28] onto [0, 55]   
    if(ieta < 0)
      iEta = ieta + 28;
    else if(ieta > 0)
      iEta = ieta + 27;
    return iEta;
  }

 int convertGenEta(double inputEta) {
   const double tpgEtaValues[27] = {
     0.087,      
     0.174, // HB and inner HE bins are 0.348 wide
     0.261,
     0.348,
     0.522,
     0.609,
     0.696,
     0.783,
     0.870,
     0.957,
     1.044,
     1.131,
     1.218,
     1.305,
     1.392,
     1.479,
     1.566,
     1.653,
     1.74,
     1.848,
     1.956, // Last two HE bins are 0.432 and 0.828 wide
     2.064,
     2.172,
     2.379,
     2.586,
     2.793,
     3
     //IGNORING HF
     //3.250, // HF bins are 0.5 wide
     //3.750,
     //4.250,
     //4.750
   };

    for (int n=1; n<29; n++){
      //std::cout<<"inputEta "<<inputEta<< " n "<< n <<" tpgEtaValues[n-1] "<< tpgEtaValues[n-1] << " abs(inputEta)<tpgEtaValues[n-1]"<<std::endl;
      if (std::fabs(inputEta)<tpgEtaValues[n-1]) {
        if(inputEta>0){ return n + 28; }
        else{ return n; }
        break;
      }
    }
    std::cout<<"OUT OF BOUNDS!!!!  inputeta: "<<inputEta<<std::endl;
    return -9;
  }

  //-pi < phi <= +pi,
  int convertGenPhi(double inputPhi){
    double posPhi[36];
    for(int n = 0; n < 36; n++)
      posPhi[n] = (0.087) * n + 0.0435;
    double negPhi[36];
    for(int n = 0; n < 36; n++)
      negPhi[n] = -3.14159 + 0.087 * n - 0.0435;

    //1 to 36 is 0 to pi
    if( 3.1416 > inputPhi && inputPhi >= 0){

      for(int n = 1; n < 36; n++){
	//std::cout<<"inputPhi "<<inputPhi<< " posPhi[n-1] "<< posPhi[n-1] << " n "<<n<<std::endl;
	if(inputPhi <= posPhi[n-1]){
	  int tpgPhi = n;
	  return tpgPhi;
	}
      }
    }

    //37 to 72 is -pi to 0
    else if(-3.1416 < inputPhi && inputPhi < 0){
      for(int n = 1; n < 36; n++)
	if(inputPhi < negPhi[n-1]){
	  int tpgPhi = n + 36;
	  return tpgPhi;
	}
    }
    std::cout<<"OUT OF BOUNDS!!!!  inputphi: "<<inputPhi<<std::endl;
    return -9;
  }

  float convertRCTEta(uint32_t inputEta) {
    const double regionEtaValues[22] = {
      -4.75,
      -4.25,
      -3.75,
      -3.25,
      -2.5,
      -1.93,
      -1.566,
      -1.218,
      -0.87,
      -0.522,
      -0.174,
      0.174,
      0.522,
      0.87,
      1.218,
      1.566,
      1.93,
      2.5,
      3.25,
      3.75,
      4.25,
      4.75
    };
    return regionEtaValues[inputEta];
  };  
  
  float convertRCTPhi(uint32_t inputPhi) {
    const double regionPhiValues[20] = {
      0.000,
      0.349,
      0.698,
      1.047,
      1.396,
      1.744,
      2.093,
      2.442,
      2.791,
      -3.14159,
      -2.791,
      -2.442,
      -2.093,
      -1.744,
      -1.396,
      -1.047,
      -0.698,
      -0.349
    };
    return regionPhiValues[inputPhi];
  };
  
  float towerEtaMap[28]=   {
    0.0435,
    0.1305, 0.2175, 0.3045, 0.3915, 0.4785,
    0.5655, 0.6525, 0.7395, 0.8265, 0.9135,
    1.0005, 1.0875, 1.1745, 1.2615, 1.3485,
    1.4355, 1.5225, 1.6095, 1.6965, 1.7835,
    1.8705, 1.9575, 2.0445, 2.217, 2.391,
    2.565, //2.739,
    2.913,
    };
  
  float towerPhiMap[72]=
    {-0.131, -0.044, 0.044, 0.131, 0.218, 0.305, 0.393, 0.480, 0.567, 0.654, 0.742, 0.829, 0.916, 1.004, 1.091, 1.178, 1.265, 1.353, 1.440, 1.527, 1.614, 1.702, 1.789, 1.876, 1.963, 2.051, 2.138, 2.225, 2.313, 2.400, 2.487, 2.574, 2.662, 2.749, 2.836, 2.923, 3.011, 3.098,
      -3.098, -3.011, -2.923, -2.836, -2.749, -2.662, -2.574, -2.487, -2.400, -2.313, -2.225, -2.138, -2.051, -1.963, -1.876, -1.789, -1.702, -1.614, -1.527, -1.440, -1.353, -1.265, -1.178, -1.091, -1.004, -0.916, -0.829, -0.742, -0.654, -0.567, -0.480, -0.393, -0.305, -0.218};

  float getRecoPhi(int iphi){
    return towerPhiMap[iphi-1];
  };

  float getRecoEta(int ieta, short zside){
    float eta = -999;
    if(ieta<0 || ieta>(28*2)){
      std::cout<<"Error!!! towereta out of bounds in triggerGeometryTools.h "<<std::endl;
      std::cout<<"ieta "<<ieta<<std::endl;
      exit(0);
    }
    
    if(zside == 1)
      eta = towerEtaMap[ieta];
    else if(zside == -1)
      eta = towerEtaMap[ieta];
    else{
      std::cout<<"Error!!! zside out of bounds in triggerGeometryTools.h "<<std::endl;
      std::cout<<"zside "<<zside<<std::endl;
      exit(0);
    }
    return eta;
  };
  
  //map of tpg eta values by ieta
  double tpgEtaMap[57] = {0};
  //map of tpg phi values by iphi
  double tpgPhiMap[73] = {0};
  //fill gen eta map and gen phi map
  void initializeEtaPhiMaps(){
    for(int i = 0; i <57; i++){
      tpgEtaMap[i] = convertGenEta(i);
    }
    for(int i = 0; i <73; i++){
      tpgPhiMap[i] = convertGenPhi(i);
    }
  }

};



#endif
