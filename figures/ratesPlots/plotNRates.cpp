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
                TString outputDir
                )
{
  assert((hists.size() == labels.size()) && (hists.size() == colors.size()));

  setTDRStyle();
  TCanvas* Tcan = new TCanvas("Tcan","", 100, 20, 1000, 800);
  TLegend* leg = new TLegend(0.55,0.15,0.90,0.45);
  applySmallerLegStyle(leg);

  Tcan->SetGrid();

  TLatex *latex = new TLatex(); 
  latex->SetNDC();
  latex->SetTextFont(42);
  latex->SetTextColor(kBlack);

  Tcan->cd();     /* Set current canvas */
  Tcan->SetFillColor(0);

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
  latex->DrawLatex(0.91, 0.960, "#scale[0.8]{0 PU}"); 

  if (!(outputName.Contains("genEta"))) {  // genPt: put legend below the efficiecy curve
    float commentaryXpos = 0.54;
    latex->DrawLatex(commentaryXpos, 0.550, "#scale[0.8]{EG Barrel}");
    latex->DrawLatex(commentaryXpos, 0.490, "#scale[0.8]{RelVal ElectronGun Pt 2 to 100}");
    latex->DrawLatex(commentaryXpos, 0.430, "#scale[0.8]{L1 p_{T} > 25, |#eta^{Gen}| < 1.4841}");
  //latex->DrawLatex(commentaryXpos, 0.700, "#scale[0.7]{}");
  //latex->DrawLatex(commentaryXpos, 0.660, bonusDescriptor);
  }  
  else { // genEta: put legend above the efficiency curve
    float commentaryXpos = 0.54;
    latex->DrawLatex(commentaryXpos, 0.900, "#scale[0.8]{EG Barrel}");
    latex->DrawLatex(commentaryXpos, 0.840, "#scale[0.8]{RelVal ElectronGun Pt 2 to 100}");
    latex->DrawLatex(commentaryXpos, 0.780, "#scale[0.8]{L1 p_{T} > 25, |#eta^{Gen}| < 1.4841}");
  } 
  Tcan->Update();


  Tcan->cd();
  Tcan->SaveAs(outputDir+outputName+".pdf");
  Tcan->SaveAs(outputDir+outputName+".png");
}
             

/*******************************************************************/

#endif
