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

//  float etaValues[95] = {-5.2665, -5.1155, -4.92125, -4.71475, -4.53875, -4.36375, -4.1895, -4.014, -3.83875, -3.664, -3.489, -3.314, -3.045, -2.958, -2.871, -2.784, -2.697, -2.61, -2.523, -2.436, -2.349, -2.262, -2.175, -2.088, -2.001, -1.914, -1.827, -1.74, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218, -1.131, -1.044, -0.957, -0.87, -0.783, -0.696, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, -0.087, 0, 0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.87, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479, 1.566, 1.653, 1.74, 1.827, 1.914, 2.001, 2.088, 2.175, 2.262, 2.349, 2.436, 2.523, 2.61, 2.697, 2.784, 2.871, 2.958, 3.045, 3.314, 3.489, 3.664, 3.83875, 4.014, 4.1895, 4.36375, 4.53875, 4.71475, 4.92125, 5.1155, 5.2665};
//  float etaValues[59] = {-2.523, -2.436, -2.349, -2.262, -2.175, -2.088, -2.001, -1.914, -1.827, -1.74, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218, -1.131, -1.044, -0.957, -0.87, -0.783, -0.696, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, -0.087, 0, 0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.87, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479, 1.566, 1.653, 1.74, 1.827, 1.914, 2.001, 2.088, 2.175, 2.262, 2.349, 2.436, 2.523 };
  float etaValues[35] = {-1.479, -1.392, -1.305, -1.218, -1.131, -1.044, -0.957, -0.87, -0.783, -0.696, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, -0.087, 0, 0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.87, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479};
  float phiValues[73] =
    {-3.142, -3.054, -2.967, -2.880, -2.793, -2.705, -2.618, -2.531, -2.443, -2.356, -2.269, -2.182, -2.094, -2.007, -1.920, -1.833, -1.745, -1.658, -1.571, -1.484, -1.396, -1.309, -1.222, -1.134, -1.047, -0.960, -0.873, -0.785, -0.698, -0.611, -0.524, -0.436, -0.349, -0.262, -0.175, -0.087,
     0.000, 0.087, 0.175, 0.262, 0.349, 0.436, 0.524, 0.611, 0.698, 0.785, 0.873, 0.960, 1.047, 1.134, 1.222, 1.309, 1.396, 1.484, 1.571, 1.658, 1.745, 1.833, 1.920, 2.007, 2.094, 2.182, 2.269, 2.356, 2.443, 2.531, 2.618, 2.705, 2.793, 2.880, 2.967, 3.054, 3.142};
  
  //eta lines
  for(int i = 1; i < 34; i++){
    TLine * line = new TLine(etaValues[i], -3.142, etaValues[i], 3.142); 
    line->SetLineColor(kGray);
    line->SetLineStyle(1);
    line->SetLineWidth(1);
    TowerLines.push_back(line);
  }

  //phi lines
  for(int i = 1; i < 72; i++){
    TLine * line = new TLine(-1.479, phiValues[i], 1.479, phiValues[i]); 
    line->SetLineColor(kGray);
    line->SetLineStyle(1);
    line->SetLineWidth(1);
    TowerLines.push_back(line);
  }

  for(size_t j = 0; j < TowerLines.size(); j++){
    TowerLines.at(j)->Draw();
  }
}

void plotEventDisplayPhaseIICaloJets(){
  
  gStyle->SetOptStat(0);
  
  TFile *f = TFile::Open("/afs/cern.ch/work/p/pdas/emulator_phase2/correlator/CMSSW_15_1_0/src/analyzer.root", "READ");

  if (!f) { return; }

  TTreeReader myReader("l1NtupleProducer/displayTree", f);
  TTreeReaderValue<vector<TLorentzVector>> vEcalTpgs(myReader, "ecalTPGs");
  TTreeReaderValue<vector<TLorentzVector>> vHcalTpgs(myReader, "hcalTPGs");
  TTreeReaderValue<vector<TLorentzVector>> vEGClusters(myReader, "egClusters");
  TTreeReaderValue<vector<TLorentzVector>> vPFClusters(myReader, "pfClusters");
  TTreeReaderValue<vector<TLorentzVector>> vOfflineJets(myReader, "offlineJets");
  TTreeReaderValue<vector<TLorentzVector>> vGctCaloJets(myReader, "gctCaloJets");
  TTreeReaderValue<vector<TLorentzVector>> vGenJets(myReader, "genJets");
  TTreeReaderValue<vector<TLorentzVector>> vGenElectrons(myReader, "genEles");
  TTreeReaderValue<int> vEvent(myReader, "event");

  //float etaValues[95] = {-5.2665, -5.1155, -4.92125, -4.71475, -4.53875, -4.36375, -4.1895, -4.014, -3.83875, -3.664, -3.489, -3.314, -3.045, -2.958, -2.871, -2.784, -2.697, -2.61, -2.523, -2.436, -2.349, -2.262, -2.175, -2.088, -2.001, -1.914, -1.827, -1.74, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218, -1.131, -1.044, -0.957, -0.87, -0.783, -0.696, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, -0.087, 0, 0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.87, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479, 1.566, 1.653, 1.74, 1.827, 1.914, 2.001, 2.088, 2.175, 2.262, 2.349, 2.436, 2.523, 2.61, 2.697, 2.784, 2.871, 2.958, 3.045, 3.314, 3.489, 3.664, 3.83875, 4.014, 4.1895, 4.36375, 4.53875, 4.71475, 4.92125, 5.1155, 5.2665};
  //float etaValues[59] = {-2.523, -2.436, -2.349, -2.262, -2.175, -2.088, -2.001, -1.914, -1.827, -1.74, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218, -1.131, -1.044, -0.957, -0.87, -0.783, -0.696, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, -0.087, 0, 0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.87, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479, 1.566, 1.653, 1.74, 1.827, 1.914, 2.001, 2.088, 2.175, 2.262, 2.349, 2.436, 2.523};
  float etaValues[35] = {-1.479, -1.392, -1.305, -1.218, -1.131, -1.044, -0.957, -0.87, -0.783, -0.696, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, -0.087, 0, 0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.87, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479};

  float phiValues[73] =
    {-3.142, -3.054, -2.967, -2.880, -2.793, -2.705, -2.618, -2.531, -2.443, -2.356, -2.269, -2.182, -2.094, -2.007, -1.920, -1.833, -1.745, -1.658, -1.571, -1.484, -1.396, -1.309, -1.222, -1.134, -1.047, -0.960, -0.873, -0.785, -0.698, -0.611, -0.524, -0.436, -0.349, -0.262, -0.175, -0.087,
     0.000, 0.087, 0.175, 0.262, 0.349, 0.436, 0.524, 0.611, 0.698, 0.785, 0.873, 0.960, 1.047, 1.134, 1.222, 1.309, 1.396, 1.484, 1.571, 1.658, 1.745, 1.833, 1.920, 2.007, 2.094, 2.182, 2.269, 2.356, 2.443, 2.531, 2.618, 2.705, 2.793, 2.880, 2.967, 3.054, 3.142};


  while (myReader.Next()) {

  // Create a new canvas
  //TCanvas *c1 = new TCanvas("c1","eta vs phi",200,10,1250,800);
  TCanvas *c1 = new TCanvas("c1","eta vs phi",200,10,950,800);
  c1->SetFillColor(0);
  c1->GetFrame()->SetFillColor(0);
  c1->GetFrame()->SetBorderSize(6);
  c1->GetFrame()->SetBorderMode(-1);
  TPad *p1 = new TPad("pad1", "pad1", 0., 0., 0.85, 1.);
  p1->Draw();
  p1->cd();

  const Int_t kUPDATE = 1000;

  // Create histograms
  TH1F   *h                = new TH1F("h","This is the eta distribution",100,-4,4);
  TH2F   *h2EcalTpgs       = new TH2F("h2EcalTpgs", "Event Display", 34, etaValues, 72, phiValues);
  TH2F   *h2HcalTpgs       = new TH2F("h2HcalTpgs", "Event Display", 34, etaValues, 72, phiValues);
  TH2F   *h2EGClusters     = new TH2F("h2EGClusters", "Event Display", 34, etaValues, 72, phiValues);
  TH2F   *h2PFClusters     = new TH2F("h2PFClusters", "Event Display", 34, etaValues, 72, phiValues);
  TH2F   *h2OfflineJets    = new TH2F("h2OfflineJets", "Event Display", 34, etaValues, 72, phiValues);
  TH2F   *h2GctCaloJets    = new TH2F("h2GctCaloJets", "Event Display", 34, etaValues, 72, phiValues);
  TH2F   *h2GenJets        = new TH2F("h2GenJets", "Event Display", 58, etaValues, 72, phiValues);
  TH2F   *h2GenElectrons   = new TH2F("h2GenElectrons", "Event Display", 34, etaValues, 72, phiValues);
  
  h->SetFillColor(48);
  int event = *vEvent;

  // Get the event number
  char name[30];
  sprintf(name,"Event %u",event);
  std::cout<<event<<std::endl;
  std::cout<<name<<std::endl;

  int ci;
  TColor *color;
  //["#1845fb", "#ff5e02", "#c91f16", "#c849a9", "#adad7d", "#86c8dd", "#578dff", "#656364"]

  // Get HCAL TPGs
  ci = TColor::GetColor("#b9ac70");
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
  h2HcalTpgs->SetFillColorAlpha(ci, 0.8);
  h2HcalTpgs->SetLineColorAlpha(ci, 0.8);
  h2HcalTpgs->GetXaxis()->SetTitle("#eta");
  h2HcalTpgs->GetYaxis()->SetTitle("#phi");
  h2HcalTpgs->SetTitle("");
  h2HcalTpgs->Draw("BOX");
  h2HcalTpgs2->SetLineColor(ci);
  h2HcalTpgs2->SetLineWidth(1);
  h2HcalTpgs2->Draw("SAME BOXL");

  //DrawCardLines();
  //DrawRegionLines();
  DrawTowerLines();
  gPad->RedrawAxis();

  // Get ECAL TPGs
  ci = TColor::GetColor("#92dadd");
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
  h2EcalTpgs->SetFillColorAlpha(ci, 0.8);
  h2EcalTpgs->SetLineColor(ci);
  h2EcalTpgs->Draw("SAME BOX");
  h2EcalTpgs2->SetLineColor(ci);
  h2EcalTpgs2->SetLineWidth(1);
  h2EcalTpgs2->Draw("SAME BOXL");

  // Get the PF clusters
  ci = TColor::GetColor("#3f90da");
  double pfMinPt = 5.0;
  if(pfMinPt > 0.){
    std::cout << "[INFO:] plotEventDisplayPhaseIICaloJets.C: do not show PF clusters with energy under "
              << pfMinPt << " GeV" << std::endl;
  }

  for (UInt_t j = 0; j < vPFClusters->size(); ++j) {
    if(vPFClusters->at(j).Pt() > pfMinPt){
      float ceta = vPFClusters->at(j).Eta();
      float cphi = vPFClusters->at(j).Phi();
      float cpt  = vPFClusters->at(j).Pt();
      h2PFClusters->Fill(ceta, cphi, cpt);
      TBox *box = new TBox(ceta-0.1305,cphi-0.1305,ceta+0.1305,cphi+0.1305);
      box->SetFillStyle(0);
      box->SetLineWidth(2);
      box->SetLineColor(ci);
      box->Draw("SAME");
      if(cpt > 10.){
        std::cout<<"vPFClusters->at(j).Pt() "<< cpt
                 <<" eta "<< ceta
                 <<" phi "<< cphi <<std::endl;
      }
    }
  }

  //h2PFClusters->SetFillStyle(1001);
  //h2PFClusters->SetFillColor(ci);
  h2PFClusters->SetLineColor(ci);
  h2PFClusters->SetLineWidth(2);
  //h2PFClusters->Draw("SAME BOXL");

  // Get the EG clusters
  ci = TColor::GetColor("#e76300");
  double egMinPt = 5.;
  if(egMinPt > 0.){
    std::cout << "[INFO:] plotEventDisplayPhaseIICaloJets.C: do not show EG clusters with energy under "
              << egMinPt << " GeV" << std::endl;
  }

  for (UInt_t j = 0; j < vEGClusters->size(); ++j) {
    if(vEGClusters->at(j).Pt() > egMinPt){
      float ceta = vEGClusters->at(j).Eta();
      float cphi = vEGClusters->at(j).Phi();
      float cpt  = vEGClusters->at(j).Pt();
      h2EGClusters->Fill(ceta, cphi, cpt);
      if(cpt > 10.){
        std::cout<<"vEGClusters->at(j).Pt() "<< cpt
                 <<" eta "<< ceta
                 <<" phi "<< cphi <<std::endl;
      }
    }
  }

  TH2F* h2EGClusters2 = (TH2F*)h2EGClusters->Clone();
  h2EGClusters->SetFillStyle(3444);
  h2EGClusters->SetFillColor(ci);
  h2EGClusters->SetLineColor(ci);
  h2EGClusters->Draw("SAME BOX");
  h2EGClusters2->SetLineColor(ci);
  h2EGClusters2->SetLineWidth(1);
  h2EGClusters2->Draw("SAME BOXL");

  // Get the offline jets
  double recoJetMinPt = 15.;
  double recoJetMaxEta = 3.0;
  std::cout << "[INFO:] plotEventDisplayPhaseIICaloJets.C: do not show offline jets with energy under "
            << recoJetMinPt << " GeV and eta greater than " << recoJetMaxEta <<std::endl;
  for (UInt_t j = 0; j < vOfflineJets->size(); ++j) {
    float ceta = vOfflineJets->at(j).Eta();
    float cphi = vOfflineJets->at(j).Phi();
    float cpt  = vOfflineJets->at(j).Pt();
    if(cpt > recoJetMinPt && fabs(ceta) < recoJetMaxEta) {
      h2OfflineJets->Fill(ceta, cphi, cpt);
      TEllipse *circ = new TEllipse(ceta,cphi,.4,.4);
      circ->SetFillStyle(0);
      circ->SetLineColor(kViolet+2);
      //circ->Draw("SAME");

      std::ostringstream strs;
      strs << cpt;
      std::string text = strs.str();
      TPaveText *tempText = new TPaveText(ceta, cphi, ceta-0.25, cphi+0.25);
      tempText->AddText(text.c_str());
      tempText->SetFillColor(0);
      tempText->SetLineColor(0);
      tempText->SetShadowColor(0);
      tempText->SetTextColor(kViolet+2);
      //tempText->Draw("SAME");
    }
  }

  h2OfflineJets->SetLineColor(kViolet+2);
  h2OfflineJets->SetLineWidth(2);
  //h2OfflineJets->Draw("SAME BOXL");

  // Get the gen jets
  double genJetMinPt = 10.;
  ci = TColor::GetColor("#717581");
  std::cout << "[INFO:] plotEventDisplayPhaseIICaloJets.C: do not show gen jets with energy under "
            << genJetMinPt << " GeV " <<std::endl;
  for (UInt_t j = 0; j < vGenJets->size(); ++j) {
    float ceta = vGenJets->at(j).Eta();
    float cphi = vGenJets->at(j).Phi();
    float cpt  = vGenJets->at(j).Pt();
    //if(cpt > genJetMinPt) {
      std::cout<<"vGenJets->at(j).Pt() "<< cpt
               <<" eta "<< ceta 
               <<" phi "<< cphi <<std::endl;
      h2GenJets->Fill(ceta, cphi, cpt);
      TEllipse *circ = new TEllipse(ceta,cphi,.4,.4);
      circ->SetFillStyle(0);
      circ->SetLineColor(ci);
      circ->SetLineWidth(2);
      if(abs(ceta) < 1.479) circ->Draw("SAME");

      std::ostringstream strs;
      strs << cpt;
      std::string text = strs.str();
      TPaveText *tempText = new TPaveText(ceta, cphi, ceta-0.25, cphi+0.25);
      tempText->AddText(text.c_str());
      tempText->SetFillColor(0);
      tempText->SetLineColor(0);
      tempText->SetShadowColor(0);
      tempText->SetTextColor(kViolet+2);
      //tempText->Draw("SAME");
    //}
  }
  //h2GenJets->SetLineColor(kViolet+2);
  h2GenJets->SetLineColor(ci);
  h2GenJets->SetLineWidth(2); 

  // Get the gen electrons
  ci = TColor::GetColor("#832db6");
  for (UInt_t j = 0; j < vGenElectrons->size(); ++j) {
    float ceta = vGenElectrons->at(j).Eta();
    float cphi = vGenElectrons->at(j).Phi();
    float cpt  = vGenElectrons->at(j).Pt();
    h2GenElectrons->Fill(ceta, cphi, cpt);
    TEllipse *circ = new TEllipse(ceta,cphi,.4,.4);
    circ->SetFillStyle(0);
    circ->SetLineColor(ci);
    circ->SetLineWidth(2);
    if(abs(ceta) < 1.479) circ->Draw("SAME");

    std::ostringstream strs;
    strs << cpt;
    std::string text = strs.str();
    TPaveText *tempText = new TPaveText(ceta, cphi, ceta-0.25, cphi+0.25);
    tempText->AddText(text.c_str());
    tempText->SetFillColor(0);
    tempText->SetLineColor(0);
    tempText->SetShadowColor(0);
    tempText->SetTextColor(kAzure+3);
    //tempText->Draw("SAME");
  }
  h2GenElectrons->SetLineColor(ci);
  h2GenElectrons->SetLineWidth(2);

  // Get the GCT jets
  ci = TColor::GetColor("#bd1f01");
  for (UInt_t j = 0; j < vGctCaloJets->size(); ++j) {
    float ceta = vGctCaloJets->at(j).Eta();
    float cphi = vGctCaloJets->at(j).Phi();
    float cpt  = vGctCaloJets->at(j).Pt();
    h2GctCaloJets->Fill(ceta, cphi, cpt);
    TBox *box = new TBox(ceta-0.3915,cphi-0.3915,ceta+0.3915,cphi+0.3915);
    box->SetFillStyle(0);
    box->SetLineWidth(2);
    box->SetLineColor(ci);
    if(abs(ceta) < 1.479) box->Draw("SAME");

    std::ostringstream strs;
    strs << cpt;
    std::string text = strs.str();
    TPaveText *tempText = new TPaveText(ceta, cphi, ceta+0.25, cphi+0.25);
    tempText->AddText(text.c_str());
    tempText->SetFillColor(0);
    tempText->SetLineColor(0);
    tempText->SetShadowColor(0);
    tempText->SetTextColor(kRed);
    //tempText->Draw("SAME");
  }

  h2GctCaloJets->SetLineColor(ci);
  h2GctCaloJets->SetLineWidth(2);
  //h2GctCaloJets->Draw("SAME BOXL");

  c1->Update();
  c1->cd();
  //gPad->Update();
  //gPad->RedrawAxis();
  float xR=0.70;
  //TLegend *l = new TLegend(xR,0.80,xR+0.30,1.0);
  TLegend *l = new TLegend(0.78,0.60,0.99,0.9);
  l->SetBorderSize(0);
  l->SetTextSize(0.03);

  TLatex *t2a = new TLatex(0.125,0.905,"#bf{CMS}");
  t2a->SetNDC();
  t2a->SetTextFont(42);
  t2a->SetTextSize(0.04);
  t2a->SetTextAlign(20);
  t2a->Draw("same");

  TLatex *t2b = new TLatex(0.34,0.9,"#bf{#it{Phase-2 Simulation Preliminary}}");
  t2b->SetNDC();
  t2b->SetTextFont(42);
  t2b->SetTextSize(0.03);
  t2b->SetTextAlign(20);
  t2b->Draw("same");

  TLatex *t2c = new TLatex(0.67,0.9,"#bf{PU 200 (14 TeV)}");
  t2c->SetNDC();
  t2c->SetTextFont(42);
  t2c->SetTextSize(0.032);
  t2c->SetTextAlign(20);
  t2c->Draw("same");

  l->AddEntry(h2EcalTpgs,      "ECAL deposit",   "F");
  l->AddEntry(h2HcalTpgs,      "HCAL deposit",     "F");
  l->AddEntry(h2EGClusters,    "EG clusters",     "F");
  l->AddEntry(h2PFClusters,    "PF clusters",     "F");
//  l->AddEntry(h2GenElectrons,  "GEN electrons",   "F");
  l->AddEntry(h2GctCaloJets,   "GCT jets",        "F");
//  l->AddEntry(h2OfflineJets,   "Offline jets",    "F");
  l->AddEntry(h2GenJets,       "Gen jets",        "F");
  l->Draw();
 
  char* saveFile = new char[200];
   
  sprintf(saveFile,"/eos/user/p/pdas/www/emulator_phase2/15_1_0/Event-%u-phase2emulator.png",event);
  c1->SaveAs(saveFile);

  sprintf(saveFile,"/eos/user/p/pdas/www/emulator_phase2/15_1_0/Event-%u-phase2emulator.pdf",event);
  c1->SaveAs(saveFile);
  }

  f->Close();
  delete f;
}
