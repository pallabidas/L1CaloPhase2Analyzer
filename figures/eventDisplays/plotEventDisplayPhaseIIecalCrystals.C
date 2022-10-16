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
#include <TLorentzVector.h>
#include <TStyle.h>
#include "TLegend.h"
#include "TEllipse.h"
#include "TPaveText.h"
#include "TLine.h"
#include <sstream>

#ifdef __MAKECINT__
#pragma link C++ class vector<float>+;
#endif

#ifdef __MAKECINT__
#pragma link C++ class vector<TLorentzVector>+;
#endif

void DrawCardLines(){
  std::vector<TLine*> cardLines;
  // float etaValues[17] = { -3, -2.107, -1.74, -1.392, -1.044, -0.696, -0.348, 0,
  //   0.348, 0.696, 1.044, 1.392, 1.74, 2.107, 3 };//0.3508
  // float phiValues[18] =
  // {-2.965, -2.617, -2.268, -1.919, -1.570, -1.221, -0.872, -0.523, -0.174, 
  //     0.174, 0.523, 0.872, 1.221, 1.570, 1.919, 2.268, 2.617, 2.965};

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
  float etaValues[13] = { -1.479, // -1.392,
                          -1.305, // -1.218, -1.131,
                          -1.044, // -0.957, -0.87,
                          -0.783, // -0.696, -0.609, 
                          -0.522, // -0.435, -0.348,
                          -0.261, // -0.174, -0.087, 
                          0,      // 0.087, 0.174,
                          0.261,  // 0.348, 0.435,
                          0.522,  // 0.609, 0.696,
                          0.783,  // 0.87, 0.957,
                          1.044,  // 1.131, 1.218,
                          1.305,  // 1.392, 
                          1.479
  };

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
  float etaValues[59] = { -3, -2.826, -2.652, -2.478, -2.304, -2.107, -2.0207, -1.93, -1.83, -1.74, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218, -1.131, -1.044, -0.957, -0.87, -0.783, -0.696, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, -0.087, 0, 0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.87, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479, 1.566, 1.653, 1.74, 1.83, 1.93, 2.0207, 2.107, 2.304, 2.478, 2.652, 2.826, 3};

  float phiValues[73] =
    {-3.142, -3.054, -2.967, -2.880, -2.793, -2.705, -2.618, -2.531, -2.443, -2.356, -2.269, -2.182, -2.094, -2.007, -1.920, -1.833, -1.745, -1.658, -1.571, -1.484, -1.396, -1.309, -1.222, -1.134, -1.047, -0.960, -0.873, -0.785, -0.698, -0.611, -0.524, -0.436, -0.349, -0.262, -0.175, -0.087,
     0.000, 0.087, 0.175, 0.262, 0.349, 0.436, 0.524, 0.611, 0.698, 0.785, 0.873, 0.960, 1.047, 1.134, 1.222, 1.309, 1.396, 1.484, 1.571, 1.658, 1.745, 1.833, 1.920, 2.007, 2.094, 2.182, 2.269, 2.356, 2.443, 2.531, 2.618, 2.705, 2.793, 2.880, 2.967, 3.054, 3.142};
  
  //eta lines
  for(int i = 0; i < 59; i++){
    TLine * line = new TLine(etaValues[i], -3.2, etaValues[i], 3.2); 
    line->SetLineColor(kGray);
    line->SetLineStyle(3);
    TowerLines.push_back(line);
  }

  //phi lines
  for(int i = 0; i < 73; i++){
    TLine * line = new TLine(-3, phiValues[i], 3, phiValues[i]); 
    line->SetLineColor(kGray);
    line->SetLineStyle(3);
    TowerLines.push_back(line);
  }

  for(size_t j = 0; j < TowerLines.size(); j++){
    TowerLines.at(j)->Draw();
  }
}


void plotEventDisplayPhaseIIecalCrystals(int iEvent){
  
  gStyle->SetOptStat(0);
  
  // float half_tower_offset = 0.04365;
  float half_tower_offset = 0.0;

  TFile *f = TFile::Open("L1EventDisplay.root", "READ");

  // Declare the center of the plot
  float etaCenter = 0.27;
  float phiCenter = -1.7889;

  if (!f) { return; }

  TTree *t = (TTree*) f->Get("l1NtupleProducer/efficiencyTree");

  std::vector<TLorentzVector> *vEcalTpgs       = 0;
  std::vector<TLorentzVector> *vHcalTpgs       = 0;
  std::vector<TLorentzVector> *vClusters       = 0;
  std::vector<TLorentzVector> *vTowers         = 0;
  std::vector<TLorentzVector> *vPFclusters     = 0;

  int event =0;

  // Create a new canvas.
  TCanvas *c1 = new TCanvas("c1","eta vs phi",200,10,700,700);
  c1->SetFillColor(0);
  c1->GetFrame()->SetFillColor(0);
  c1->GetFrame()->SetBorderSize(6);
  c1->GetFrame()->SetBorderMode(-1);

  const Int_t kUPDATE = 1000;

  TBranch *bEvent            = 0;      
  TBranch *bEcalTpgs         = 0;
  TBranch *bHcalTpgs         = 0;
  TBranch *bClusters = 0;
  TBranch *bTowers   = 0;
  TBranch *bPFclusters  = 0;

  t->SetBranchAddress("event",&event,&bEvent);

  t->SetBranchAddress("ecalTPGs",&vEcalTpgs,&bEcalTpgs);
  t->SetBranchAddress("hcalTPGs",&vHcalTpgs,&bHcalTpgs);
  t->SetBranchAddress("ecalClusters",&vClusters,&bClusters);
  t->SetBranchAddress("caloTowers",&vTowers,&bTowers);
  t->SetBranchAddress("caloPFClusters",&vPFclusters,&bPFclusters);

  // Create one histograms
  TH1F   *h                = new TH1F("h","This is the eta distribution",100,-4,4);
  TH2F   *h2EcalTpgs       = new TH2F("h2EcalTpgs","Event Display",(34*5), //(90*2), //64*2
				      -1.4841, 1.4841,
				      (72*5), // (144*2),
				      -3.142,3.142);
  TH2F   *h2HcalTpgs       = new TH2F("h2HcalTpgs","Event Display",34,
				      -1.4841, 1.4841,
				      72,
				      -3.142,3.142);
  TH2F   *h2L1Clusters  = new TH2F("h2L1Clusters","Event Display", (34*5), //(90*2),
				   -1.4841, 1.4841,
				   (72*5),//(144*2),
				   -3.142,3.142); 
  TH2F   *h2PFclusters  = new TH2F("h2PFclusters", "Event Display", 34,
                                   -1.4841, 1.4841,
                                   72,
                                   -3.142,3.142);
  TH2F   *h2L1Towers    = new TH2F("h2L1Towers", "Event Display", 34,
				   -1.4841, 1.4841,
				   72,
				   -3.142,3.142);
  
  h->SetFillColor(48);

  int i = iEvent;
  Long64_t tentry = t->LoadTree(i);
  std::cout<<"i "<<i<< " tentry "<< tentry << std::endl;
  bEvent->GetEntry(tentry);
  bEcalTpgs->GetEntry(tentry);
  bHcalTpgs->GetEntry(tentry);
  bClusters->GetEntry(tentry);
  bTowers->GetEntry(tentry);
  bPFclusters->GetEntry(tentry);

  //get the event number
  char name[30];
  sprintf(name,"Event %u",event);
  std::cout<<event<<std::endl;
  std::cout<<name<<std::endl;

  // Get ECAL TPGs
  double ecalMinPt = 0.;
  if(ecalMinPt > 0.){
    std::cout << "[INFO:] plotEventDisplayPhaseIIecalCrystals.C: do not show ECAL TPGs with energy under "
              << ecalMinPt << " GeV" << std::endl;
  }

  for (UInt_t j = 0; j < vEcalTpgs->size(); ++j) {
    if(vEcalTpgs->at(j).Pt() > ecalMinPt) {
      float ceta = vEcalTpgs->at(j).Eta();
      float cphi = vEcalTpgs->at(j).Phi();
      float cpt  = vEcalTpgs->at(j).Pt();

      h2EcalTpgs->Fill(ceta, cphi, cpt);
  
      std::cout<<"vEcalTpgs->at(j).Pt() "<< cpt
             <<" eta "<< ceta
             <<" phi "<< cphi <<std::endl;
    }
  }

  // Get HCAL TPGs
  double hcalMinPt = 0.;
  if(hcalMinPt > 0.){
    std::cout << "[INFO:] plotEventDisplayPhaseIIecalCrystals.C: do not show HCAL TPGs with energy under "
              << hcalMinPt << " GeV" << std::endl;
  }

  for (UInt_t j = 0; j < vHcalTpgs->size(); ++j) {
    if(vHcalTpgs->at(j).Pt() > hcalMinPt) {
      float ceta = vHcalTpgs->at(j).Eta();
      float cphi = vHcalTpgs->at(j).Phi();
      float cpt  = vHcalTpgs->at(j).Pt();

      h2HcalTpgs->Fill(ceta, cphi, cpt);

      if ((ceta > (etaCenter - 0.25)) && (ceta < (etaCenter + 0.25))
          && (cphi > (phiCenter - 0.25)) && (cphi < (phiCenter + 0.25))) {

        std::cout<<"vHcalTpgs->at(j).Pt() "<< cpt
                 <<" eta "<< ceta
                 <<" phi "<< cphi <<std::endl;
      }
    }
  }

  // Get the clusters
  double clusterMinPt = 0.;
  if(clusterMinPt > 0.){
    std::cout << "[INFO:] plotEventDisplayPhaseIIecalCrystals.C: do not show EG clusters with energy under "
              << clusterMinPt << " GeV" << std::endl;
  }

  for (UInt_t j = 0; j < vClusters->size(); ++j) {
    if(vClusters->at(j).Pt() > clusterMinPt){
      float ceta = vClusters->at(j).Eta();
      float cphi = vClusters->at(j).Phi();
      float cpt  = vClusters->at(j).Pt();

      h2L1Clusters->Fill(ceta, cphi, cpt);

      if ((ceta > (etaCenter - 0.25)) && (ceta < (etaCenter + 0.25))
          && (cphi > (phiCenter - 0.25)) && (cphi < (phiCenter + 0.25)) ) {
          std::cout<<"vClusters->at(j).Pt() "<< cpt
                   <<" eta "<< ceta
                   <<" phi "<< cphi <<std::endl;
      }
    }
  }

  // Get the towers
  double towerMinPt = 0.;
  if(towerMinPt > 0.){
    std::cout << "[INFO:] plotEventDisplayPhaseIIecalCrystals.C: do not show GCT towers with energy under "
              << towerMinPt << " GeV" << std::endl;
  }
  
  for (UInt_t j = 0; j < vTowers->size(); ++j) {
    if(vTowers->at(j).Pt() > towerMinPt){
      float ceta = vTowers->at(j).Eta();
      float cphi = vTowers->at(j).Phi();
      float cpt  = vTowers->at(j).Pt();

      h2L1Towers->Fill(ceta, cphi, cpt);

      if ((ceta > (etaCenter - 0.25)) && (ceta < (etaCenter + 0.25))
         && (cphi > (phiCenter - 0.25)) && (cphi < (phiCenter + 0.25))) {

        std::cout<<"vTowers->at(j).Pt() "<< cpt
                 <<" eta "<< ceta
                 <<" phi "<< cphi <<std::endl;
      }
    }
  }

  // Plot the HCAL TPGs (first, clone to plot the border)                                                 
  TH2F* h2HcalTpgs2 = (TH2F*)h2HcalTpgs->Clone();
  h2HcalTpgs->SetFillStyle(1001);
  h2HcalTpgs->SetFillColorAlpha(kSpring+10, 0.8);
  h2HcalTpgs->SetLineColorAlpha(kSpring+10, 0.8);
  h2HcalTpgs->GetXaxis()->SetTitle("#eta");
  h2HcalTpgs->GetYaxis()->SetTitle("#phi");
  h2HcalTpgs->SetTitle(name);

  h2HcalTpgs->GetXaxis()->SetRangeUser(etaCenter - 0.25, etaCenter + 0.25);
  h2HcalTpgs->GetYaxis()->SetRangeUser(phiCenter - 0.25, phiCenter + 0.25);

  h2HcalTpgs->Draw("BOX");
  h2HcalTpgs->Draw("SAME BOX");
  h2HcalTpgs2->SetLineColor(kSpring+10);
  h2HcalTpgs2->SetLineWidth(1);
  h2HcalTpgs2->Draw("SAME BOXL");

  DrawCardLines();
  DrawRegionLines();
  DrawTowerLines();

  // Plot the ECAL crystals (TPGs)
  TH2F* h2EcalTpgs2 = (TH2F*)h2EcalTpgs->Clone();
  h2EcalTpgs->SetFillStyle(1001);
  h2EcalTpgs->SetFillColorAlpha(kPink+1, 0.8);
  h2EcalTpgs->SetLineColorAlpha(kPink+1, 0.8);
  h2EcalTpgs->Draw("SAME BOX");
  h2EcalTpgs2->SetLineColor(kPink+1);
  h2EcalTpgs2->SetLineWidth(1);
  h2EcalTpgs2->Draw("SAME BOXL");

  // Plot the towers
  TH2F* h2L1Towers2 = (TH2F*)h2L1Towers->Clone();
  h2L1Towers->SetFillStyle(3444);
  h2L1Towers->SetFillColor(kGreen+3);
  h2L1Towers->SetLineColor(kGreen+3);
  h2L1Towers->Draw("SAME BOX");
  h2L1Towers2->SetLineColor(kGreen+3);
  h2L1Towers2->SetLineWidth(1);
  h2L1Towers2->Draw("SAME BOXL");

  // Plot the clusters
  TH2F* h2L1Clusters2 = (TH2F*)h2L1Clusters->Clone();
  h2L1Clusters->SetFillStyle(1001);
  h2L1Clusters->SetFillColorAlpha(kRed, 0.8);
  h2L1Clusters->SetLineColorAlpha(kRed, 0.8);
  h2L1Clusters->Draw("SAME BOX");
  h2L1Clusters2->SetLineColor(kRed);
  h2L1Clusters2->SetLineWidth(1);
  h2L1Clusters2->Draw("SAME BOXL");

  float xR=0.70;
  TLegend *l = new TLegend(xR,0.80,xR+0.30,1.0);
  l->SetTextSize(0.035);

  l->AddEntry(h2EcalTpgs,   "ECAL Crystals",   "F");
  l->AddEntry(h2HcalTpgs,   "HCAL Towers",     "F");
  l->AddEntry(h2L1Clusters, "EG Clusters", "F");
  l->AddEntry(h2L1Towers,   "GCT Towers",   "F");
  l->Draw();

  char* saveFile = new char[100];

  sprintf(saveFile,"/eos/user/s/skkwan/phase2RCTDevel/events/Event-%u-current_emulator_test.png",event);
  c1->SaveAs(saveFile);

  sprintf(saveFile,"/eos/user/s/skkwan/phase2RCTDevel/events/Event-%u-current_emulator_test.pdf",event);
  c1->SaveAs(saveFile);

}
