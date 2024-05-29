#include <vector>
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2.h"
#include "TGraph2D.h"
#include "TGraph.h"
#include "TBenchmark.h"
#include "TRandom.h"
#include "TSystem.h"
#include "TLorentzVector.h"
#include <TStyle.h>
#include "TLegend.h"
#include "TEllipse.h"
#include "TPaveText.h"
#include "TLine.h"
#include "TEllipse.h"
#include <sstream>
//#include "Math/VectorUtil_Cint.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

#ifdef __MAKECINT__
#pragma link C++ class std::vector<float>+;
#endif

#ifdef __MAKECINT__
#pragma link C++ class std::vector<TLorentzVector>+;
#endif

void DrawCardLines(){
  std::vector<TLine*> cardLines;

  float etaValues[3] = { -1.479, 0, 1.479 };

  float phiValues[19] =
    { -3.142, -2.793, -2.443, -2.094, -1.745, -1.396, -1.047, -0.698, -0.349, 0.000, 
      0.349, 0.698, 1.047, 1.396, 1.745, 2.094, 2.443, 2.793, 3.142};
  
  //eta lines
  for(int i = 0; i < 3; i++){
    TLine * line = new TLine(etaValues[i], -3.2, etaValues[i], 3.2); 
    line->SetLineColor(kRed);
    line->SetLineStyle(1);
    line->SetLineWidth(2);
    cardLines.push_back(line);
  }

  //phi lines
  for(int i = 0; i < 19; i++){
    TLine * line = new TLine(-3, phiValues[i], 3, phiValues[i]); 
    line->SetLineColor(kRed);
    line->SetLineStyle(1);
    line->SetLineWidth(2);
    cardLines.push_back(line);
  }

  for(size_t j = 0; j < cardLines.size(); j++){
    cardLines.at(j)->Draw();
  }
}


/*
 * Draw ECAL region lines.
 */
void DrawRegionLines(){

  std::vector<TLine*> RegionLines;
  float etaValues[13] = { -1.479, -1.305, -1.044, -0.783, -0.522, -0.261, 0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479};

  float phiValues[19] =
    { -3.142, -2.793, -2.443, -2.094, -1.745, -1.396, -1.047, -0.698, -0.349, 0.000,
      0.349, 0.698, 1.047, 1.396, 1.745, 2.094, 2.443, 2.793, 3.142};

  //eta lines
  for(int i = 0; i < 13; i++){
    TLine * line = new TLine(etaValues[i], -3.2, etaValues[i], 3.2); 
    line->SetLineColor(kBlue);
    line->SetLineStyle(1);
    RegionLines.push_back(line);
  }

  //phi lines
  for(int i = 0; i < 19; i++){
    TLine * line = new TLine(-3, phiValues[i], 3, phiValues[i]); 
    line->SetLineColor(kBlue);
    line->SetLineStyle(1);
    RegionLines.push_back(line);
  }

  for(size_t j = 0; j < RegionLines.size(); j++){
    RegionLines.at(j)->Draw();
  }


}

/*
 * Draw tower lines.
 */
void DrawTowerLines(){
  std::vector<TLine*> TowerLines;

  float etaValues[95] = {-5.2665, -5.1155, -4.92125, -4.71475, -4.53875, -4.36375, -4.1895, -4.014, -3.83875, -3.664, -3.489, -3.314, -3.045, -2.958, -2.871, -2.784, -2.697, -2.61, -2.523, -2.436, -2.349, -2.262, -2.175, -2.088, -2.001, -1.914, -1.827, -1.74, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218, -1.131, -1.044, -0.957, -0.87, -0.783, -0.696, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, -0.087, 0, 0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.87, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479, 1.566, 1.653, 1.74, 1.827, 1.914, 2.001, 2.088, 2.175, 2.262, 2.349, 2.436, 2.523, 2.61, 2.697, 2.784, 2.871, 2.958, 3.045, 3.314, 3.489, 3.664, 3.83875, 4.014, 4.1895, 4.36375, 4.53875, 4.71475, 4.92125, 5.1155, 5.2665};

  float phiValues[73] =
    {-3.142, -3.054, -2.967, -2.880, -2.793, -2.705, -2.618, -2.531, -2.443, -2.356, -2.269, -2.182, -2.094, -2.007, -1.920, -1.833, -1.745, -1.658, -1.571, -1.484, -1.396, -1.309, -1.222, -1.134, -1.047, -0.960, -0.873, -0.785, -0.698, -0.611, -0.524, -0.436, -0.349, -0.262, -0.175, -0.087,
     0.000, 0.087, 0.175, 0.262, 0.349, 0.436, 0.524, 0.611, 0.698, 0.785, 0.873, 0.960, 1.047, 1.134, 1.222, 1.309, 1.396, 1.484, 1.571, 1.658, 1.745, 1.833, 1.920, 2.007, 2.094, 2.182, 2.269, 2.356, 2.443, 2.531, 2.618, 2.705, 2.793, 2.880, 2.967, 3.054, 3.142};
  
  //eta lines
  for(int i = 1; i < 94; i++){
    TLine * line = new TLine(etaValues[i], -3.142, etaValues[i], 3.142); 
    line->SetLineColor(kGray);
    line->SetLineStyle(1);
    line->SetLineWidth(1);
    TowerLines.push_back(line);
  }

  //phi lines
  for(int i = 1; i < 72; i++){
    TLine * line = new TLine(-5.2665, phiValues[i], 5.2665, phiValues[i]); 
    line->SetLineColor(kGray);
    line->SetLineStyle(1);
    line->SetLineWidth(1);
    TowerLines.push_back(line);
  }

  for(size_t j = 0; j < TowerLines.size(); j++){
    TowerLines.at(j)->Draw();
  }
}

void plotEventDisplayPhaseIIPFclusters(){
  
  gStyle->SetOptStat(0);
  
  TFile *f = TFile::Open("/afs/cern.ch/work/p/pdas/emulator_phase2/calojet/sep2023/njets/CMSSW_14_1_0_pre1/src/L1Trigger/L1CaloPhase2Analyzer/test/analyzer.root", "READ");

  if (!f) { return; }

  TTreeReader myReader("l1NtupleProducer/displayTree", f);
  TTreeReaderValue<vector<TLorentzVector>> vEcalTpgs(myReader, "ecalTPGs");
  TTreeReaderValue<vector<TLorentzVector>> vHcalTpgs(myReader, "hcalTPGs");
  TTreeReaderValue<vector<TLorentzVector>> vTowers(myReader, "gctTowers");
  TTreeReaderValue<vector<TLorentzVector>> vHgcalTowers(myReader, "hgcalTowers");
  TTreeReaderValue<vector<TLorentzVector>> vHfTowers(myReader, "hfTowers");
  TTreeReaderValue<vector<TLorentzVector>> vPFclusters(myReader, "caloPFClusters");
//  TTreeReaderValue<vector<TLorentzVector>> vOfflineJets(myReader, "offlineJets");
//  TTreeReaderValue<vector<TLorentzVector>> vGctCaloJets(myReader, "gctCaloJets");
//  TTreeReaderValue<vector<TLorentzVector>> vGenJets(myReader, "genJets");
  TTreeReaderValue<int> vEvent(myReader, "event");

  float etaValues[95] = {-5.2665, -5.1155, -4.92125, -4.71475, -4.53875, -4.36375, -4.1895, -4.014, -3.83875, -3.664, -3.489, -3.314, -3.045, -2.958, -2.871, -2.784, -2.697, -2.61, -2.523, -2.436, -2.349, -2.262, -2.175, -2.088, -2.001, -1.914, -1.827, -1.74, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218, -1.131, -1.044, -0.957, -0.87, -0.783, -0.696, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, -0.087, 0, 0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.87, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479, 1.566, 1.653, 1.74, 1.827, 1.914, 2.001, 2.088, 2.175, 2.262, 2.349, 2.436, 2.523, 2.61, 2.697, 2.784, 2.871, 2.958, 3.045, 3.314, 3.489, 3.664, 3.83875, 4.014, 4.1895, 4.36375, 4.53875, 4.71475, 4.92125, 5.1155, 5.2665};

  float phiValues[73] =
    {-3.142, -3.054, -2.967, -2.880, -2.793, -2.705, -2.618, -2.531, -2.443, -2.356, -2.269, -2.182, -2.094, -2.007, -1.920, -1.833, -1.745, -1.658, -1.571, -1.484, -1.396, -1.309, -1.222, -1.134, -1.047, -0.960, -0.873, -0.785, -0.698, -0.611, -0.524, -0.436, -0.349, -0.262, -0.175, -0.087,
     0.000, 0.087, 0.175, 0.262, 0.349, 0.436, 0.524, 0.611, 0.698, 0.785, 0.873, 0.960, 1.047, 1.134, 1.222, 1.309, 1.396, 1.484, 1.571, 1.658, 1.745, 1.833, 1.920, 2.007, 2.094, 2.182, 2.269, 2.356, 2.443, 2.531, 2.618, 2.705, 2.793, 2.880, 2.967, 3.054, 3.142};


  while (myReader.Next()) {

  // Create a new canvas
  TCanvas *c1 = new TCanvas("c1","eta vs phi",200,10,1250,800);
  c1->SetFillColor(0);
  c1->GetFrame()->SetFillColor(0);
  c1->GetFrame()->SetBorderSize(6);
  c1->GetFrame()->SetBorderMode(-1);
  TPad *p1 = new TPad("pad1", "pad1", 0., 0., 0.9, 1.);
  p1->Draw();
  p1->cd();

  const Int_t kUPDATE = 1000;

  // Create histograms
  TH1F   *h                = new TH1F("h","This is the eta distribution",100,-4,4);
  TH2F   *h2EcalTpgs       = new TH2F("h2EcalTpgs", "Event Display", 94, etaValues, 72, phiValues);
  TH2F   *h2HcalTpgs       = new TH2F("h2HcalTpgs", "Event Display", 94, etaValues, 72, phiValues);
  TH2F   *h2HgcalTowers    = new TH2F("h2HgcalTowers", "Event Display", 94, etaValues, 72, phiValues);
  TH2F   *h2HfTowers       = new TH2F("h2HfTowers", "Event Display", 94, etaValues, 72, phiValues);
  TH2F   *h2L1Towers       = new TH2F("h2L1Towers", "Event Display", 94, etaValues, 72, phiValues);
  TH2F   *h2PFclusters     = new TH2F("h2PFclusters", "Event Display", 94, etaValues, 72, phiValues);
//  TH2F   *h2OfflineJets    = new TH2F("h2OfflineJets", "Event Display", 94, etaValues, 72, phiValues);
//  TH2F   *h2GctCaloJets    = new TH2F("h2GctCaloJets", "Event Display", 94, etaValues, 72, phiValues);
//  TH2F   *h2GenJets        = new TH2F("h2GenJets", "Event Display", 94, etaValues, 72, phiValues);
  
  h->SetFillColor(48);
  int event = *vEvent;

  // Get the event number
  char name[30];
  sprintf(name,"Event %u",event);
  std::cout<<event<<std::endl;
  std::cout<<name<<std::endl;

  // Get HCAL TPGs
  double hcalMinPt = 0.5;
  if(hcalMinPt > 0.){
    std::cout << "[INFO:] plotEventDisplayPhaseIICaloJets.C: do not show HCAL TPGs with energy under "
              << hcalMinPt << " GeV" << std::endl;
  }

  for (UInt_t j = 0; j < vHcalTpgs->size(); ++j) {
    if(vHcalTpgs->at(j).Pt() > hcalMinPt) {
      float ceta = vHcalTpgs->at(j).Eta();
      float cphi = vHcalTpgs->at(j).Phi();
      float cpt  = vHcalTpgs->at(j).Pt();

      h2HcalTpgs->Fill(ceta, cphi, cpt);

      if(cpt > 10.){

        std::cout<<"vHcalTpgs->at(j).Pt() "<< cpt
                 <<" eta "<< ceta
                 <<" phi "<< cphi <<std::endl;
      }
    }
  }

  TH2F* h2HcalTpgs2 = (TH2F*)h2HcalTpgs->Clone();
  h2HcalTpgs->SetFillStyle(1001);
  h2HcalTpgs->SetFillColorAlpha(kSpring+10, 0.8);
  h2HcalTpgs->SetLineColorAlpha(kSpring+10, 0.8);
  h2HcalTpgs->GetXaxis()->SetTitle("#eta");
  h2HcalTpgs->GetYaxis()->SetTitle("#phi");
  h2HcalTpgs->SetTitle("");
  h2HcalTpgs->Draw("BOX");
  h2HcalTpgs2->SetLineColor(kSpring+10);
  h2HcalTpgs2->SetLineWidth(1);
  h2HcalTpgs2->Draw("SAME BOXL");

  //DrawCardLines();
  //DrawRegionLines();
  DrawTowerLines();
  gPad->RedrawAxis();

  // Get ECAL TPGs
  double ecalMinPt = 0.5;
  if(ecalMinPt > 0.){
    std::cout << "[INFO:] plotEventDisplayPhaseIICaloJets.C: do not show ECAL TPGs with energy under "
              << ecalMinPt << " GeV" << std::endl;
  }

  for (UInt_t j = 0; j < vEcalTpgs->size(); ++j) {
    if(vEcalTpgs->at(j).Pt() > ecalMinPt) {
      float ceta = vEcalTpgs->at(j).Eta();
      float cphi = vEcalTpgs->at(j).Phi();
      float cpt  = vEcalTpgs->at(j).Pt();

      h2EcalTpgs->Fill(ceta, cphi, cpt);

      if(cpt > 10.){
        std::cout<<"vEcalTpgs->at(j).Pt() "<< cpt
               <<" eta "<< ceta
               <<" phi "<< cphi <<std::endl;
      }
    }
  }

  TH2F* h2EcalTpgs2 = (TH2F*)h2EcalTpgs->Clone();
  h2EcalTpgs->SetFillStyle(1001);
  h2EcalTpgs->SetFillColorAlpha(kPink+1, 0.8);
  h2EcalTpgs->SetLineColorAlpha(kPink+1, 0.8);
  h2EcalTpgs->Draw("SAME BOX");
  h2EcalTpgs2->SetLineColor(kPink+1);
  h2EcalTpgs2->SetLineWidth(1);
  h2EcalTpgs2->Draw("SAME BOXL");

  // Get the GCTintTowers
  double towerMinPt = 0.5;
  if(towerMinPt > 0.){
    std::cout << "[INFO:] plotEventDisplayPhaseIICaloJets.C: do not show GCT towers with energy under "
              << towerMinPt << " GeV" << std::endl;
  }

  for (UInt_t j = 0; j < vTowers->size(); ++j) {
    if(vTowers->at(j).Pt() > towerMinPt){
      float ceta = vTowers->at(j).Eta();
      float cphi = vTowers->at(j).Phi();
      float cpt  = vTowers->at(j).Pt();

      h2L1Towers->Fill(ceta, cphi, cpt);

      if(cpt > 10.){
        std::cout<<"vTowers->at(j).Pt() "<< cpt
                 <<" eta "<< ceta
                 <<" phi "<< cphi <<std::endl;
      }
    }
  }

  TH2F* h2L1Towers2 = (TH2F*)h2L1Towers->Clone();
  h2L1Towers->SetFillStyle(3444);
  h2L1Towers->SetFillColor(kAzure+10);
  h2L1Towers->SetLineColor(kAzure+10);
  h2L1Towers->Draw("SAME BOX");
  h2L1Towers2->SetLineColor(kAzure+10);
  h2L1Towers2->SetLineWidth(1);
  h2L1Towers2->Draw("SAME BOXL");

  // Get the HGCAL towers
  double hgcalMinPt = 1.0;
  if(hgcalMinPt > 0.){
    std::cout << "[INFO:] plotEventDisplayPhaseIICaloJets.C: do not show HGCAL towers with energy under "
              << hgcalMinPt << " GeV" << std::endl;
  }

  for (UInt_t j = 0; j < vHgcalTowers->size(); ++j) {
    if(vHgcalTowers->at(j).Pt() > hgcalMinPt){
      float ceta = vHgcalTowers->at(j).Eta();
      float cphi = vHgcalTowers->at(j).Phi();
      float cpt  = vHgcalTowers->at(j).Pt();
      h2HgcalTowers->Fill(ceta, cphi, cpt);

      if(cpt > 10.){
        std::cout<<"vHgcalTowers->at(j).Pt() "<< cpt
                 <<" eta "<< ceta
                 <<" phi "<< cphi <<std::endl;
      }
    }
  }

  h2HgcalTowers->SetFillStyle(1001);
  h2HgcalTowers->SetFillColor(kBlue-9);
  h2HgcalTowers->SetLineColor(kBlue-9);
  h2HgcalTowers->Draw("SAME BOX");
  h2HgcalTowers->SetLineColor(kBlue-9);
  h2HgcalTowers->SetLineWidth(1);
  h2HgcalTowers->Draw("SAME BOXL");

  // Get the HF towers
  double hfMinPt = 1.0;
  if(hfMinPt > 0.){
    std::cout << "[INFO:] plotEventDisplayPhaseIICaloJets.C: do not show HF towers with energy under "
              << hfMinPt << " GeV" << std::endl;
  }

  for (UInt_t j = 0; j < vHfTowers->size(); ++j) {
    if(vHfTowers->at(j).Pt() > hfMinPt){
      float ceta = vHfTowers->at(j).Eta();
      float cphi = vHfTowers->at(j).Phi();
      float cpt  = vHfTowers->at(j).Pt();
      h2HfTowers->Fill(ceta, cphi, cpt);

      if(cpt > 10.){
        std::cout<<"vHfTowers->at(j).Pt() "<< cpt
                 <<" eta "<< ceta
                 <<" phi "<< cphi <<std::endl;
      }
    }
  }

  h2HfTowers->SetFillStyle(1001);
  h2HfTowers->SetFillColor(kBlue-6);
  h2HfTowers->SetLineColor(kBlue-6);
  h2HfTowers->Draw("SAME BOX");
  h2HfTowers->SetLineColor(kBlue-6);
  h2HfTowers->SetLineWidth(1);
  h2HfTowers->Draw("SAME BOXL");

  // Get the PF clusters
  double pfclusterMinPt = 0.;
  if(pfclusterMinPt > 0.){
    std::cout << "[INFO:] plotEventDisplayPhaseIIPFclusters.C: do not show GCT towers with energy under "
              << pfclusterMinPt << " GeV" << std::endl;
  }

  for (UInt_t j = 0; j < vPFclusters->size(); ++j) {
    if(vPFclusters->at(j).Pt() > pfclusterMinPt){
      float ceta = vPFclusters->at(j).Eta();
      float cphi = vPFclusters->at(j).Phi();
      float cpt  = vPFclusters->at(j).Pt();
      h2PFclusters->Fill(ceta, cphi, cpt);

//      if(cpt > 10.){
//        std::cout<<"vPFclusters->at(j).Pt() "<< cpt
//                 <<" eta "<< ceta
//                 <<" phi "<< cphi <<std::endl;
//      }
    }
  }

  h2PFclusters->SetLineColor(kViolet-6);
  h2PFclusters->SetLineWidth(2);
  h2PFclusters->Draw("SAME BOXL");

//  // Get the Offline Jets
//  double recoJetMinPt = 15.;
//  double recoJetMaxEta = 3.0;
//  std::cout << "[INFO:] plotEventDisplayPhaseIICaloJets.C: do not show offline jets with energy under "
//            << recoJetMinPt << " GeV and eta greater than " << recoJetMaxEta <<std::endl;
//  for (UInt_t j = 0; j < vOfflineJets->size(); ++j) {
//    float ceta = vOfflineJets->at(j).Eta();
//    float cphi = vOfflineJets->at(j).Phi();
//    float cpt  = vOfflineJets->at(j).Pt();
//    if(cpt > recoJetMinPt && fabs(ceta) < recoJetMaxEta) {
//      h2OfflineJets->Fill(ceta, cphi, cpt);
//      TEllipse *circ = new TEllipse(ceta,cphi,.4,.4);
//      circ->SetFillStyle(0);
//      circ->SetLineColor(kViolet+2);
//      //circ->Draw("SAME");
//
//      std::ostringstream strs;
//      strs << cpt;
//      std::string text = strs.str();
//      TPaveText *tempText = new TPaveText(ceta, cphi, ceta-0.25, cphi+0.25);
//      tempText->AddText(text.c_str());
//      tempText->SetFillColor(0);
//      tempText->SetLineColor(0);
//      tempText->SetShadowColor(0);
//      tempText->SetTextColor(kViolet+2);
//      //tempText->Draw("SAME");
//    }
//  }
//
//  h2OfflineJets->SetLineColor(kViolet+2);
//  h2OfflineJets->SetLineWidth(2);
//  //h2OfflineJets->Draw("SAME BOXL");
//
//  double genJetMinPt = 10.;
//  std::cout << "[INFO:] plotEventDisplayPhaseIICaloJets.C: do not show gen jets with energy under "
//            << genJetMinPt << " GeV " <<std::endl;
//  for (UInt_t j = 0; j < vGenJets->size(); ++j) {
//    float ceta = vGenJets->at(j).Eta();
//    float cphi = vGenJets->at(j).Phi();
//    float cpt  = vGenJets->at(j).Pt();
//    if(cpt > genJetMinPt) {
//      std::cout<<"vGenJets->at(j).Pt() "<< cpt
//               <<" eta "<< ceta 
//               <<" phi "<< cphi <<std::endl;
//      h2GenJets->Fill(ceta, cphi, cpt);
//      TEllipse *circ = new TEllipse(ceta,cphi,.4,.4);
//      circ->SetFillStyle(0);
//      circ->SetLineColor(kViolet+2);
//      circ->Draw("SAME");
//
//      std::ostringstream strs;
//      strs << cpt;
//      std::string text = strs.str();
//      TPaveText *tempText = new TPaveText(ceta, cphi, ceta-0.25, cphi+0.25);
//      tempText->AddText(text.c_str());
//      tempText->SetFillColor(0);
//      tempText->SetLineColor(0);
//      tempText->SetShadowColor(0);
//      tempText->SetTextColor(kViolet+2);
//      tempText->Draw("SAME");
//    }
//  }
//  h2GenJets->SetLineColor(kViolet+2);
//  h2GenJets->SetLineWidth(2); 
//
//  // Get the Calo Jets
//  for (UInt_t j = 0; j < vGctCaloJets->size(); ++j) {
//    float ceta = vGctCaloJets->at(j).Eta();
//    float cphi = vGctCaloJets->at(j).Phi();
//    float cpt  = vGctCaloJets->at(j).Pt();
//    h2GctCaloJets->Fill(ceta, cphi, cpt);
//    TEllipse *circ = new TEllipse(ceta,cphi,.4,.4);
//    circ->SetFillStyle(0);
//    circ->SetLineColor(kRed);
//    //circ->Draw("SAME");
//
//    std::ostringstream strs;
//    strs << cpt;
//    std::string text = strs.str();
//    TPaveText *tempText = new TPaveText(ceta, cphi, ceta+0.25, cphi+0.25);
//    tempText->AddText(text.c_str());
//    tempText->SetFillColor(0);
//    tempText->SetLineColor(0);
//    tempText->SetShadowColor(0);
//    tempText->SetTextColor(kRed);
//    tempText->Draw("SAME");
//  }
//
//  h2GctCaloJets->SetLineColor(kRed);
//  h2GctCaloJets->SetLineWidth(2);
//  h2GctCaloJets->Draw("SAME BOXL");

  c1->Update();
  c1->cd();
  //gPad->Update();
  //gPad->RedrawAxis();
  float xR=0.70;
  //TLegend *l = new TLegend(xR,0.80,xR+0.30,1.0);
  TLegend *l = new TLegend(0.82,0.60,0.99,0.9);
  l->SetTextSize(0.03);

//  TLatex *t2a = new TLatex(0.46,0.9," #bf{CMS} #it{Phase-2 Simulation Preliminary}                14 TeV, 200PU   ");
//  t2a->SetNDC();
//  t2a->SetTextFont(42);
//  t2a->SetTextSize(0.045);
//  t2a->SetTextAlign(20);
//  t2a->Draw("same");

  TLatex *t2a = new TLatex(0.125,0.905,"#bf{CMS}");
  t2a->SetNDC();
  t2a->SetTextFont(42);
  t2a->SetTextSize(0.045);
  t2a->SetTextAlign(20);
  t2a->Draw("same");

  TLatex *t2b = new TLatex(0.3,0.9,"#bf{Phase-2 Simulation Preliminary}");
  t2b->SetNDC();
  t2b->SetTextFont(42);
  t2b->SetTextSize(0.032);
  t2b->SetTextAlign(20);
  t2b->Draw("same");

  TLatex *t2c = new TLatex(0.73,0.9,"#bf{14 TeV, 200PU}");
  t2c->SetNDC();
  t2c->SetTextFont(42);
  t2c->SetTextSize(0.04);
  t2c->SetTextAlign(20);
  t2c->Draw("same");

  l->AddEntry(h2EcalTpgs,    "ECAL Crystals",   "F");
  l->AddEntry(h2HcalTpgs,    "HCAL Towers",     "F");
  l->AddEntry(h2L1Towers,    "GCT Towers",      "F");
  l->AddEntry(h2HgcalTowers, "HGCal Towers",    "F");
  l->AddEntry(h2HfTowers,    "HF Towers",       "F");
  l->AddEntry(h2PFclusters,  "PF clusters",     "F");
  //l->AddEntry(h2GctCaloJets, "GCT Jets",    "F");
  //l->AddEntry(h2OfflineJets, "Offline Jets",    "F");
  //l->AddEntry(h2GenJets,     "Gen Jets",        "F");
  l->Draw();
 
  char* saveFile = new char[200];
   
  sprintf(saveFile,"/eos/user/p/pdas/www/emulator_phase2/14_1_0/pfcluster_test/Event-%u-phase2emulator.png",event);
  c1->SaveAs(saveFile);

  sprintf(saveFile,"/eos/user/p/pdas/www/emulator_phase2/14_1_0/pfcluster_test/Event-%u-phase2emulator.pdf",event);
  c1->SaveAs(saveFile);
  }

  f->Close();
  delete f;
}
