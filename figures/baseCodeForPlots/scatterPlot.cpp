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

#include "singleDistribution.C"

#ifndef SCATTER_PLOT_C_INCL
#define SCATTER_PLOT_C_INCL
 

 
/* Generate a single scatter plot. N.B.: has some more argument than singleDistribution.C, namely "yLabel"
   and the y-axes properties.
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
int scatterPlot(TString name, TString variable, TString cut, TString legend,
                            TString treePath, TString inputDirectory, TString outputDirectory,
                            TString xLabel, TString yLabel, 
                            TString bonusDescriptor,
                            int xbins, float xlow, float xhigh,
                            int ybins, float ylow, float yhigh){ 
 
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
 
  TH2F *hist = new TH2F("hist","hist", xbins, xlow, xhigh, ybins, ylow, yhigh);
  tree->Draw(variable+">>+hist", cut);

  hist->SetMarkerColor(kBlack);
  hist->SetMarkerStyle(kFullDotLarge);
  hist->SetLineWidth(2);
  hist->SetLineColor(kBlack);


  hist->Draw("HIST"); 

  
  hist->GetXaxis()->SetTitle(xLabel);
  hist->GetXaxis()->SetLabelSize(0.04);
  hist->GetXaxis()->SetNdivisions(-505);


  hist->GetYaxis()->SetTitle(yLabel);
  hist->GetYaxis()->SetLabelSize(0.04);
  hist->GetYaxis()->SetNdivisions(505, kTRUE);


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


  float commentaryXpos = 0.23;
  latex->DrawLatex(0.80, 0.920, "#scale[0.8]{0 PU}");
  latex->DrawLatex(commentaryXpos, 0.840, "#scale[0.6]{EG Barrel}");
  latex->DrawLatex(commentaryXpos, 0.800, "#scale[0.6]{RelVal ElectronGun Pt 2 to 100}");
  latex->DrawLatex(commentaryXpos, 0.760, "#scale[0.6]{v0 ECAL propagation (no ECAL stitching),}");
  latex->DrawLatex(commentaryXpos, 0.720, "#scale[0.6]{|#eta^{Gen}| < 1.4841}");
  latex->DrawLatex(commentaryXpos, 0.660, bonusDescriptor);

  Tcan->Update();

  Tcan->SaveAs(outputDirectory+name+".png");
  Tcan->SaveAs(outputDirectory+name+".pdf");
 
  delete Tcan;

  return 1;

}



#endif
