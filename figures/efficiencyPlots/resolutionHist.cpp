/*******************************************************************/
/* efficiencyHist.cpp                                              */
/* Author: Stephanie Kwan                                          */
/*******************************************************************/

#include "TH1F.h"
#include "TH1.h"
#include "TH1F.h"

#include "TAxis.h"
#include "TChain.h"

#include "THStack.h"
#include "TLegend.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TTree.h"

#include <assert.h>
#include <cmath>
#include <iostream>
#include <math.h>
#include <string>
#include <TMath.h>
#include <vector>

#include "../baseCodeForPlots/tdrstyle.C"
#include "../baseCodeForPlots/comparisonPlots.C"

#ifndef RESOLUTION_HIST_CPP_INCL
#define RESOLUTION_HIST_CPP_INCL

/*******************************************************************/

void plotNResolutions(std::vector<TH1F*> graphs, 
             std::vector<TString> labels,
             std::vector<int> colors,
             TString xAxisLabel,
             TString legendName,
             TString outputName,
             TString outputDir
             )
{
  assert((graphs.size() == labels.size()) && (graphs.size() == colors.size()));

  setTDRStyle();
  TCanvas* Tcan = new TCanvas("Tcan","", 100, 20, 1000, 800);
  TLegend* leg = new TLegend(0.65,0.15,0.90,0.45);
  applySmallerLegStyle(leg);

  Tcan->SetGrid();

  TLatex *latex = new TLatex(); 
  latex->SetNDC();
  latex->SetTextFont(42);
  latex->SetTextColor(kBlack);

  Tcan->cd();     /* Set current canvas */
  Tcan->SetFillColor(0);


  std::vector<TH1F*>::iterator itGraph;
  std::vector<TString>::iterator itLabel;
  std::vector<int>::iterator itColor;
  
  TH1F* histDummy;
  for (itGraph = graphs.begin(), itLabel = labels.begin(), itColor = colors.begin();
       itGraph != graphs.end();
       itGraph++, itLabel++, itColor++)
    {
      if (itGraph == graphs.begin()) // only do this once 
   {
     histDummy = new TH1F(**itGraph);
   }
      
      // De-reference the iterator to get the TH1F*
      (*itGraph)->SetMarkerColor(*itColor);
      (*itGraph)->SetMarkerStyle(kFullCircle);
      (*itGraph)->SetLineWidth(2);
      (*itGraph)->SetLineColor(*itColor);
      (*itGraph)->SetMarkerSize(2);
    }

  histDummy->SetMarkerColor(0);
  histDummy->SetLineColor(0);

  histDummy->Draw("");

  for (itGraph = graphs.begin(); itGraph != graphs.end(); itGraph++)
    {
      (*itGraph)->Draw("hist same");
    }

  //histDummy->GetXaxis()->SetTitle("Resolution vs Gen Electron p_{T} [GeV]");
  histDummy->GetXaxis()->SetTitle("Resolution vs Gen Pion p_{T} [GeV]");
  histDummy->GetYaxis()->SetTitle("#entries");
  histDummy->GetXaxis()->SetTitleSize(0.06); // default is 0.03                                                                    
  /* Set y-axis limits */  
  //histDummy->GetYaxis()->SetRangeUser(0.0, 1.5);
  // histDummy->GetYaxis()->SetRangeUser(0.8, 1.02);

  /* Customize legend */
  for (itGraph = graphs.begin(), itLabel = labels.begin();
       itGraph != graphs.end();
       itGraph++, itLabel++)
    {
      leg->AddEntry(*itGraph, *itLabel,  "l");
    }
  leg->Draw();


  // Default to RCT label, use GCT if not
  TString emuLabel = "#scale[1.0]{#bf{CMS}} #scale[0.8]{#it{Phase 2 GCT emulator}}";  
  if (outputName.Contains("RCT")) {
    emuLabel = "#scale[1.0]{#bf{CMS}} #scale[0.8]{#it{Phase 2 RCT emulator}}";  
  }
  latex->DrawLatex(0.16, 0.960, emuLabel); 
  latex->DrawLatex(0.91, 0.960, "#scale[0.8]{0 PU}"); 

//  if (!(outputName.Contains("genEta")) && !(outputName.Contains("genPhi"))) {  // genPt: put legend below the efficiecy curve
//    float commentaryXpos = 0.54;
//    latex->DrawLatex(commentaryXpos, 0.550, "#scale[0.8]{EG Barrel}");
//    latex->DrawLatex(commentaryXpos, 0.490, "#scale[0.8]{RelVal ElectronGun Pt 2 to 100}");
//    latex->DrawLatex(commentaryXpos, 0.430, "#scale[0.8]{L1 p_{T} > 25, |#eta^{Gen}| < 1.4841}");
//  //latex->DrawLatex(commentaryXpos, 0.700, "#scale[0.7]{}");
//  //latex->DrawLatex(commentaryXpos, 0.660, bonusDescriptor);
//  }  
//  else { // genEta: put legend above the efficiency curve
    float commentaryXpos = 0.54;
    latex->DrawLatex(commentaryXpos, 0.900, "#scale[0.8]{EG Barrel}");
    //latex->DrawLatex(commentaryXpos, 0.840, "#scale[0.8]{RelVal ElectronGun Pt 2 to 100}");
    latex->DrawLatex(commentaryXpos, 0.840, "#scale[0.8]{SinglePion Pt 0 to 200}");
    latex->DrawLatex(commentaryXpos, 0.780, "#scale[0.8]{L1 p_{T} > 2, |#eta^{Gen}| < 1.4841}");
//  } 
  Tcan->Update();


  Tcan->cd();
  Tcan->SaveAs(outputDir+outputName+".pdf");
  Tcan->SaveAs(outputDir+outputName+".png");

  Tcan->Close();
  delete Tcan;
}
             

/*******************************************************************/

#endif
