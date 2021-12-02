#include "TAxis.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TFormula.h"
#include "TGaxis.h"
#include "TH1.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TLegend.h"
#include "THStack.h"
#include "TMath.h"
#include "TPad.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TTree.h"

#include <cmath>
#include <iostream>
#include <math.h>
#include <string>

#include "tdrstyle.C"
#include "CMS_lumi.h"

#ifndef COMPARISON_PLOTS_INCL
#define COMPARISON_PLOTS_INCL
 
/* Apply template style to a TPad* pad1. */
void applyPadStyle(TPad* pad1){
  pad1->SetFillColor(0);
  pad1->Draw();  pad1->cd();  pad1->SetLeftMargin(0.2);  pad1->SetBottomMargin(0.13); pad1->SetRightMargin(0.1);
  //pad1->SetGrid(); 
  pad1->SetGrid(10,10); 
}

/* Apply legend style to a TLegend *leg. */
void applyLegStyle(TLegend *leg){
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetHeader("");
  leg->SetShadowColor(0);
}
 
/* Apply smaller legend style to a TLegend *leg. */
void applySmallerLegStyle(TLegend *leg){
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.044);
  leg->SetTextFont(42);
  leg->SetHeader("");
  leg->SetShadowColor(0);
}

/* Generate a comparison plot of "variable", using the gen-level genCut, and selecting fakes using fakeCut. treePath specifies the tree in the ROOT file to use.
   The ROOT file is located at inputDirectory. The resulting plots are written to outputDirectory, with filename including "name". The histogram has (bins)
   number of bins and ranges from integers low to high. */
int comparisonPlots(TString variable, TString genCut, TString fakeCut, TString treePath, TString inputDirectory, TString outputDirectory, TString name, int bins, int low, int high){ 
 
  //gROOT->LoadMacro("CMS_lumi.C");
  //gROOT->ProcessLine(".L ~/Documents/work/Analysis/PhaseIIStudies/2018/tdrstyle.C");
  setTDRStyle();
 
  //int bins = 30;
  //int low = 0;
  //int high = 15;
 
  //  TFile* tauFile = new TFile("dummy");
  TCanvas* Tcan = new TCanvas("Tcan","", 100, 20, 800, 600);
  Tcan->cd();     /* Set current canvas */
  Tcan->SetFillColor(0);
  //TPad* pad1 = new TPad("pad1","The pad",0,0.0,0.98,1);
  //applyPadStyle(pad1);
 
  TLegend *leg = new TLegend(0.60,0.75,0.85,0.9);
  applyLegStyle(leg);
 
  TFile *file = new TFile(inputDirectory);
 
  if(!file->IsOpen()||file==0){
    std::cout<<"ERROR FILE "<< inputDirectory + name<<" NOT FOUND; EXITING"<<std::endl;
    return 0;
  }
 
  //gStyle->SetOptFit(0);
  //gStyle->SetOptStat(0);
  //gStyle->SetEndErrorSize(0);
  //gStyle->SetErrorX(0.5);
 
  TTree* tree = (TTree*)file->Get(treePath);
  if(tree == 0){
    std::cout<<"ERROR Tau Tree is "<< tree<<" NOT FOUND; EXITING"<<std::endl;
    return 0;
  }
 
  TH1F *True = new TH1F("True","True",bins,low,high);
  tree->Draw(variable+">>+True", genCut);
  TH1F *Fake = new TH1F("Fake","Fake",bins,low,high);
  tree->Draw(variable+">>+Fake", fakeCut);  
 
  /* Compute the ratio by cloning the True histogram, subtracting the Fake values,
     and dividing it by the Fake value. */
  /*
  TH1F *ratio = (TH1F*)True->Clone();
  ratio->Add(Fake,-1);
  ratio->Divide(Fake);
  */

  True->SetMarkerColor(0);
  True->SetFillStyle(1001);
  True->SetFillColorAlpha(kBlue+2, 0.1);
  True->SetLineWidth(1);
  True->SetLineColor(kBlue+2);

  Fake->SetMarkerColor(0);
  Fake->SetFillStyle(1001);
  Fake->SetFillColorAlpha(kRed+2, 0.1);
  Fake->SetLineWidth(1);
  Fake->SetLineColor(kRed+2);

  Fake->Scale(1/Fake->Integral());

  True->Scale(1/True->Integral());
  Tcan->SetLogy();

  // NoIso->Draw();
  Fake->Draw("HIST");  
  True->Draw("HIST same"); 
  
  Fake->GetXaxis()->SetTitle(variable);
  Fake->GetYaxis()->SetTitle("A.U.");

  /*
  float max = 10;
  if(Fake->GetXaxis()->GetBinCenter( Fake->GetMaximumBin() ) > True->GetXaxis()->GetBinCenter( True->GetMaximumBin() ) )
    max = Fake->GetXaxis()->GetBinCenter(Fake->GetMaximumBin());
  else
    max = True->GetXaxis()->GetBinCenter(True->GetMaximumBin());
   
    std::cout<<"max: "<<max<<std::endl;
  Fake->SetMaximum(max/60);
  */

  //leg->AddEntry(NoIso,"No Isolation","l");
  //  leg->AddEntry(True,"#tau_{h} Gen-Vis p_{T}>20 GeV","l");
  //  leg->AddEntry(Fake,"Fake Background","l");
  leg->AddEntry(True,"L1 #tau_{h}, L1 p_{T}>0 GeV","l");
  leg->AddEntry(Fake,"Fake Background","l");
  leg->Draw();
 
  Tcan->cd();

  //TPad* pad2 = new TPad("pad2","The lower pad",0,0,0.98,0.25);
  //applyPadStyle(pad2);
  //pad2->cd();
  //pad2->SetGrid(0,0); 
 
  //ratio->Draw("p");

 
  Tcan->cd();
  Tcan->SetLogy();
  // Tcan->SaveAs(outputDirectory+name+".pdf");
  Tcan->SaveAs(outputDirectory+name+".png");
 
  delete Tcan;

  return 1;

}

#endif
