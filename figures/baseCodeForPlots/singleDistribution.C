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
  //pad1->SetGrid(10,10); 
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
   name: name of the output .png/.pdf
   variable: the string that will be passed to TTree->Draw()
   cut: baseline cut
   legend: legend on the x-axis (usually a description of the variable), LaTeX allowed
   treePath: path to the TTree inside the ROOT file, e.g. folderName/efficiencyTree
   inputDirectory: relative/absolute path to the ROOT file
   outputDirectory: relative/absolute path to the folder where .pngs and .pdfs will be written
   bins: number of bins in the histo
   low: min x-range
   high: max x-range
   ymax: if non-zero positive value, histogram will be capped at this max value on the y-axis
         (start with a -99 dummy value, then adjust/fine tune when you know what it will look like)
   */
int singleDistributionPlots(TString name, TString variable, TString cut, TString legend,
                            TString treePath, TString inputDirectory, TString outputDirectory,
			                      TString xLabel, TString bonusDescriptor,
                            int bins, float low, float high, float ymax){ 
 
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
 
  TH1F *hist = new TH1F("hist","hist", bins,low, high);
  tree->Draw(variable+">>+hist", cut);

  hist->SetMarkerColor(0);
  hist->SetLineWidth(2);
  hist->SetLineColor(kBlack);

  hist->Scale(1/hist->Integral());
  //  Tcan->SetLogy();

  hist->Draw("HIST"); 

  
  hist->GetXaxis()->SetTitle(xLabel);
  hist->GetXaxis()->SetLabelSize(0.04);
  hist->GetXaxis()->SetNdivisions(-505);


  hist->GetYaxis()->SetTitle("Events_{bin}/Events_{tot}");
  hist->GetYaxis()->SetLabelSize(0.04);
  hist->GetYaxis()->SetNdivisions(505, kTRUE);

  // Features specific to this analyzer: if the plot is rct_deltaR, draw a line at x = crystalSize
  float crystalSize = 0.01746;
  if (variable.Contains("deltaR")) {

    TLine *line = new TLine(crystalSize, 0, crystalSize, 1.10* hist->GetMaximum());
    line->SetLineColor(kRed);
    line->SetLineStyle(kDashed);
    line->SetLineWidth(5);
    line->Draw();

    TLatex *l = new TLatex(); 
    l->SetNDC();
    l->SetTextFont(42);
    l->SetTextColor(kRed);
    l->DrawLatex(0.25, 0.75, "#scale[0.8]{Crystal size}");
    Tcan->Update();
  }

  // More specifics: if the plot is the pT difference %, draw a line at x = 0
  float zeroDiff = 0;
  if (name.Contains("pT_fractional_diff")) {
        TLine *line = new TLine(zeroDiff, 0, zeroDiff, 1.05* hist->GetMaximum()); // 1.10* hist->GetMaximum());
    line->SetLineColor(kBlue);
    line->SetLineStyle(kDashed);
    line->SetLineWidth(3);
    line->Draw();

    TLatex *l = new TLatex(); 
    l->SetNDC();
    l->SetTextFont(42);
    l->SetTextColor(kBlue);
    l->DrawLatex(0.538, 0.091, "#scale[0.9]{0}");
    Tcan->Update();

    // Make line at -0.10 for debugging studies
    float debugDiff = -0.1;
    TLine *lineDebug = new TLine(debugDiff, 0, debugDiff, 1.05* hist->GetMaximum()); // 1.10* hist->GetMaximum());
    lineDebug->SetLineColor(kBlack);
    lineDebug->SetLineStyle(kDashed);
    lineDebug->SetLineWidth(3);
    lineDebug->Draw();
    // Make line at -0.10 for debugging studies
    TLatex *lDebug = new TLatex(); 
    lDebug->SetNDC();
    lDebug->SetTextFont(42);
    lDebug->SetTextColor(kBlack);
    lDebug->DrawLatex(0.420, 0.50, "#scale[0.9]{-0.1}");
    Tcan->Update();
  }

  // More specifics: overlay the eta boundaries
    // More specifics: if the plot is the pT difference %, draw a line at x = 0
  if (name.Contains("Eta")) {
    std:vector<float> etaBoundaries{-1.479, -1.3096, -1.0476, -0.7857, -0.5238, -0.2619, 0, 0.2619, 0.5238, 0.7857, 1.0476, 1.3096, 1.479};
    for (float b : etaBoundaries) {

  //      TLine *line = new TLine(b, 0.012, b, 1.05* hist->GetMaximum()); // 1.10* hist->GetMaximum());
      TLine *line = new TLine(b, 0.012, b, 1.05* hist->GetMaximum()); // 1.10* hist->GetMaximum());
      line->SetLineColor(kBlue);
      line->SetLineStyle(kDashed);
      line->SetLineWidth(3);
      line->Draw();
      Tcan->Update();
    }
  }

  // Finicky ymax values
  // For specific distributions, give more room at the top
  if (name.Contains("Eta") || name.Contains("Phi") || name.Contains("pT_fractional_diff")) {  
    float max = hist->GetMaximum(); 
    hist->SetMaximum(max * 1.5);
  }
  if (ymax > 0) {
    hist->SetMaximum(ymax);
  }


  // // leg->AddEntry(hist, legend,"l");
  // // leg->Draw();
 
  Tcan->cd();
 


  TLatex *latex = new TLatex(); 
  latex->SetNDC();
  latex->SetTextFont(42);
  latex->SetTextColor(kBlack);

  // Default to RCT label, use GCT if not
  TString emuLabel = "#scale[1.0]{#bf{CMS}} #scale[0.8]{#it{Phase 2 GCT emulator}}";  
  if (variable.Contains("rct")) {
    emuLabel = "#scale[1.0]{#bf{CMS}} #scale[0.8]{#it{Phase 2 RCT emulator}}";  
  }
  latex->DrawLatex(0.19, 0.920, emuLabel);  
  // latex->DrawLatex(0.19, 0.920, "#scale[1.0]{#bf{CMS}} #scale[0.8]{#it{Phase 2 RCT emulator}}"); 


  float commentaryXpos = 0.47;
  latex->DrawLatex(0.80, 0.920, "#scale[0.8]{0 PU}");
  latex->DrawLatex(commentaryXpos, 0.820, "#scale[0.6]{EG Barrel}");
  latex->DrawLatex(commentaryXpos, 0.780, "#scale[0.6]{RelVal ElectronGun Pt 2 to 100}");
  latex->DrawLatex(commentaryXpos, 0.740, "#scale[0.6]{v0 ECAL propagation,}");

  if (bonusDescriptor != "") { // 1.4841
    latex->DrawLatex(commentaryXpos, 0.700, "#scale[0.6]{|#eta^{Gen}| < 1.4841, " + bonusDescriptor + "}");
  }
  else {
    latex->DrawLatex(commentaryXpos, 0.700, "#scale[0.6]{|#eta^{Gen}| < 1.4841}");
  }
  Tcan->Update();

  Tcan->SaveAs(outputDirectory+name+".png");
  Tcan->SaveAs(outputDirectory+name+".pdf");
 
  delete Tcan;

  return 1;

}



#endif
