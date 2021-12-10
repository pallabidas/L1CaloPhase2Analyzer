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

#ifndef SINGLE_PLOTS_INCL
#define SINGLE_PLOTS_INCL
 
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

 
/* Read in and plot an EXISTING TH1F:
    th1fName :  the name of the TH1F inside the file (usually the name of the variable).
    xLabel   :  the x-axis label on the plot.
    paveText :  text to go in a TPave box. 
    isAU     : if True, normalize the plot so the area under the curve is 1, and set y-axis title to "A.U.".
               if False, do not normalize and set y-axis title to "Counts".
    histPath :  the path to the TH1F inside the file (usually the name of the folder in the ROOT file).
    inputDirectory : the path to the ROOT file (~/myRepo/input.root).
    outputDirectory: the directory where the output plots will be saved.
    useLogy:  if True, use log y-axis.
*/

int plotTH1F(TString th1fName, 
             TString xLabel, TString paveText,
             bool isAU,
             TString histPath, 
			       TString inputDirectory,
             TString outputDirectory,
             bool useLogy){ 
 
  //gROOT->LoadMacro("CMS_lumi.C");
  //gROOT->ProcessLine(".L ~/Documents/work/Analysis/PhaseIIStudies/2018/tdrstyle.C");
//  gROOT->ProcessLine(".L ../baseCodeForPlots/singleVariablePlots.C");
  setTDRStyle();

  TCanvas* Tcan = new TCanvas("Tcan","", 100, 20, 800, 600);


  Tcan->cd();     /* Set current canvas */
  Tcan->SetFillColor(0);
  //TPad* pad1 = new TPad("pad1","The pad",0,0.0,0.98,1);
  //applyPadStyle(pad1);

  if (useLogy) gPad->SetLogy();
 
  TLegend *leg = new TLegend(0.60,0.75,0.85,0.9);
  applyLegStyle(leg);
 
  TFile *file = new TFile(inputDirectory);
 
  if(!file->IsOpen()||file==0){
    std::cout<<"ERROR FILE "<< inputDirectory<<" NOT FOUND; EXITING"<<std::endl;
    return 0;
  }
 
  //gStyle->SetOptFit(0);
  //gStyle->SetOptStat(0);
  //gStyle->SetEndErrorSize(0);
  //gStyle->SetErrorX(0.5);
 
 
  TH1F* hist = 0;
  file->GetObject(histPath+th1fName, hist);

  hist->SetMarkerColor(0);
  // hist->SetFillStyle(1001);
  //  hist->SetFillColorAlpha(kBlue+2, 0.1);
  hist->SetLineWidth(1);
  hist->SetLineColor(kBlue+2);

  ///// DRAWING

  // If normalizing the histogram, y-axis should be A.U.
  if (isAU) {
    hist->Scale(1/hist->Integral());
    hist->GetYaxis()->SetTitle("A.U.");
  }
  else {
    hist->GetYaxis()->SetTitle("Counts");
  }
  hist->GetYaxis()->SetTitleSize(0.05);
  
  
  // X-axis label
  hist->GetXaxis()->SetTitle(xLabel);
  hist->GetXaxis()->SetTitleSize(0.05);


  
  hist->Draw("HIST"); 

  // Legend
  // leg->AddEntry(hist, legLabel, "l");
  // leg->Draw();

  // Draw the TPaveText
  Tcan->cd();
  TPaveText *pt = new TPaveText(0.6, 0.85, 0.96, 0.90, "NDCNB"); 
  pt->SetFillColor(kWhite);
  pt->AddText(paveText);
  pt->SetTextColor(kBlack);
  pt->Draw();
 
  //  Save the plot
  Tcan->cd();
  Tcan->SaveAs(outputDirectory+th1fName+".png");
  Tcan->SaveAs(outputDirectory+th1fName+".pdf");
 
  delete Tcan;


  return 1;

}

#endif
