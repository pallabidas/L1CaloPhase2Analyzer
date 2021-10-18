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

#ifndef SINGLE_DISTRIBUTION_C_INCL
#define SINGLE_DISTRIBUTION_C_INCL 
 
/* Apply template style to a TPad* pad1. */
void applyPadStyle(TPad* pad1){
  pad1->SetFillColor(0);
  pad1->Draw();  pad1->cd();  pad1->SetLeftMargin(0.2);  pad1->SetBottomMargin(0.13); pad1->SetRightMargin(0.1);
  pad1->SetTopMargin(0.1);
  //pad1->SetGrid(); 
  pad1->SetGrid(10,10); 
}

/* Apply legend style to a TLegend *leg. */
void applyLegStyle(TLegend *leg){
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.06);
  leg->SetTextFont(42);
  leg->SetHeader("");
  leg->SetShadowColor(0);
}
 
/* Generate a single distribution plot.
   treePath specifies the tree in the ROOT file to use.
   The ROOT file is located at inputDirectory. The resulting plots are written to outputDirectory, with filename including "variable". The histogram has (bins)
   number of bins and ranges from integers low to high.
   "legend" is the legend label, "xLabel" is the x-axis label. */
int singleDistributionPlots(TString variable, TString cut, TString legend, TString treePath, TString inputDirectory, TString outputDirectory,
			    TString xLabel, int bins, float low, float high){ 
 
  gROOT->LoadMacro("/Users/stephaniekwan/Documents/Phase2L1Calo/phase2-l1Calo-analyzer/figures/baseCodeForPlots/CMS_lumi.C");
  //gROOT->ProcessLine(".L ~/Documents/work/Analysis/PhaseIIStudies/2018/tdrstyle.C");
  setTDRStyle();

 
  TCanvas* Tcan = new TCanvas("Tcan","", 200, 40, 1600, 1600);
  Tcan->cd();     /* Set current canvas */
  Tcan->SetFillColor(0);
  TPad* pad1 = new TPad("pad1","The pad",0,0.0,0.98,1);
  applyPadStyle(pad1);

 
  // TLegend *leg = new TLegend(0.20,0.85,0.9,0.95);
  // applyLegStyle(leg);
 
  TFile *file = new TFile(inputDirectory);
 
  if(!file->IsOpen()||file==0){
    std::cout<<"ERROR FILE "<< inputDirectory <<" NOT FOUND; EXITING"<<std::endl;
    return 0;
  }
 
  // //gStyle->SetOptFit(0);
  // //gStyle->SetOptStat(0);
  // //gStyle->SetEndErrorSize(0);
  // //gStyle->SetErrorX(0.5);
 
  TTree* tree = (TTree*)file->Get(treePath);
  if(tree == 0){
    std::cout<<"ERROR: " << treePath <<" NOT FOUND; EXITING"<<std::endl;
    return 0;
  }
 
  TH1F *hist = new TH1F("hist","hist", bins,low,high);
  tree->Draw(variable+">>+hist", cut);

  hist->SetMarkerColor(0);
  hist->SetLineWidth(4);
  hist->SetLineColor(kBlack);

  hist->Scale(1/hist->Integral());
  //  Tcan->SetLogy();

  hist->Draw("HIST"); 

  
  hist->GetXaxis()->SetTitle(xLabel);
  hist->GetXaxis()->SetLabelSize(0.04);


  hist->GetYaxis()->SetTitle("Events_{bin}/Events_{tot}");
  hist->GetYaxis()->SetLabelSize(0.04);

  // Features specific to this analyzer: if the plot is deltaR, draw a line at the y = crystalSize
  float crystalSize = 0.01746;
  if (variable == "deltaR") {

    TLine *line = new TLine(crystalSize, 0, crystalSize, hist->GetMaximum());
    line->SetLineColor(kRed);
    line->SetLineStyle(kDashed);
    line->SetLineWidth(5);
    line->Draw();


  }

  // //  leg->AddEntry(hist,"#tau_{h} Gen-Vis p_{T}>20 GeV","l");


  // // leg->AddEntry(hist, legend,"l");
  // // leg->Draw();
 
  Tcan->cd();
 

  TLatex *latex = new TLatex(); 
  latex->SetNDC();
  latex->SetTextFont(42);
  latex->SetTextColor(kBlack);
  latex->DrawLatex(0.19, 0.920, "#scale[1.2]{#bf{CMS}} #it{Preliminary}");
  Tcan->Update();

  Tcan->SaveAs(outputDirectory+variable+".png");
  Tcan->SaveAs(outputDirectory+variable+".pdf");
 
  delete Tcan;

  return 1;

}

#endif
