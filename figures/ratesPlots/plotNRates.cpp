/*******************************************************************/
/* plotRates.cpp                                                   */
/* Author: Stephanie Kwan                                          */
/*******************************************************************/

#include "TGraphAsymmErrors.h"
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

#ifndef PLOTRATES_CPP_INCL
#define PLOTRATES_CPP_INCL

/*******************************************************************/

void plotNRates(std::vector<TH1F*> hists, 
                std::vector<TString> labels,
                std::vector<int> colors,
                float xMin, float xMax,
                float yMin, float yMax,
                TString xAxisLabel,
                TString outputName,
                TString outputDir,
                bool useLogy
                )
{
  assert((hists.size() == labels.size()) && (hists.size() == colors.size()));

  setTDRStyle();
  TCanvas* Tcan = new TCanvas("Tcan","", 100, 20, 1000, 800);
//  TLegend* leg = new TLegend(0.55,0.15,0.90,0.45); // bottom right corner
  TLegend* leg = new TLegend(0.55,0.65,0.90,0.95);
  applySmallerLegStyle(leg);

  Tcan->SetGrid();

  TLatex *latex = new TLatex(); 
  latex->SetNDC();
  latex->SetTextFont(42);
  latex->SetTextColor(kBlack);

  Tcan->cd();     /* Set current canvas */
  Tcan->SetFillColor(0);

  if (useLogy) gPad->SetLogy();


  // Make canvas larger so y-axis label won't run off the side
  Tcan->SetLeftMargin(0.18);

  std::vector<TH1F*>::iterator itHist;
  std::vector<TString>::iterator itLabel;
  std::vector<int>::iterator itColor;
  
  TH1F* histDummy;
  for (itHist = hists.begin(), itLabel = labels.begin(), itColor = colors.begin();
       itHist != hists.end();
       itHist++, itLabel++, itColor++)
    {
      if (itHist == hists.begin()) // only do this once 
   {
     histDummy = new TH1F(**itHist);
   }
      
      // De-reference the iterator to get the TH1F*
      (*itHist)->SetLineWidth(3);
      (*itHist)->SetLineColor(*itColor);
    }


  histDummy->Draw("");

  for (itHist = hists.begin(); itHist != hists.end(); itHist++)
    {
      (*itHist)->Draw("SAME");
    }

  histDummy->GetXaxis()->SetRangeUser(xMin, xMax);
  histDummy->GetYaxis()->SetRangeUser(yMin, yMax);
  histDummy->GetXaxis()->SetTitle(xAxisLabel);
  histDummy->GetYaxis()->SetTitle("Rate [kHz]");
  histDummy->GetXaxis()->SetTitleSize(0.06); // default is 0.03                                                                    
  histDummy->GetYaxis()->SetTitleSize(0.06); // default is 0.03 

  histDummy->GetYaxis()->SetMaxDigits(6); // Suppress scientific notation on y-axis because it clashes with the header
  // histDummy->GetYaxis()->SetNoExponent(kFALSE);
  /* Customize legend */
  //  leg->SetHeader(legendName); 
  for (itHist = hists.begin(), itLabel = labels.begin();
       itHist != hists.end();
       itHist++, itLabel++)
    {
      leg->AddEntry(*itHist, *itLabel,  "l");
    }
  leg->Draw();


  // Default to RCT label, use GCT if not
  TString emuLabel = "#scale[1.0]{#bf{CMS}} #scale[0.8]{#it{Phase 2 GCT emulator}}";  
  if (outputName.Contains("RCT")) {
    emuLabel = "#scale[1.0]{#bf{CMS}} #scale[0.8]{#it{Phase 2 RCT emulator}}";  
  }
  latex->DrawLatex(0.17, 0.960, emuLabel); 
  latex->DrawLatex(0.89, 0.960, "#scale[0.8]{200 PU}"); 

  float commentaryXpos = 0.54;
  latex->DrawLatex(commentaryXpos, 0.850, "#scale[0.8]{EG Barrel, MinBias 200 PU}");
  latex->DrawLatex(commentaryXpos, 0.790, "#scale[0.8]{Phase 2 HLT TDR Winter20}");

  Tcan->Update();


  Tcan->cd();
  Tcan->SaveAs(outputDir+outputName+".pdf");
  Tcan->SaveAs(outputDir+outputName+".png");
}
             

/*******************************************************************/

#endif
