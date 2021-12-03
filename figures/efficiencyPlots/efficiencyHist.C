/*******************************************************************/
/* efficiencyHist.cpp                                              */
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

#ifndef EFFICIENCY_HIST_CPP_INCL
#define EFFICIENCY_HIST_CPP_INCL

/*******************************************************************/

void plotNEfficiencies(std::vector<TGraphAsymmErrors*> graphs, 
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
  TLegend* leg = new TLegend(0.60,0.75,0.85,0.9);
  Tcan->SetGrid();

  Tcan->cd();     /* Set current canvas */
  Tcan->SetFillColor(0);
  leg = new TLegend(0.55, 0.2, 0.90, 0.5);
  applyLegStyle(leg);

  std::vector<TGraphAsymmErrors*>::iterator itGraph;
  std::vector<TString>::iterator itLabel;
  std::vector<int>::iterator itColor;
  
  TGraphAsymmErrors* histDummy;
  for (itGraph = graphs.begin(), itLabel = labels.begin(), itColor = colors.begin();
       itGraph != graphs.end();
       itGraph++, itLabel++, itColor++)
    {
      if (itGraph == graphs.begin()) // only do this once 
   {
     histDummy = new TGraphAsymmErrors(**itGraph);
   }
      
      // De-reference the iterator to get the TGraphAsymmErrors*
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
      (*itGraph)->Draw("P");
    }

  histDummy->GetXaxis()->SetTitle(xAxisLabel);
  histDummy->GetYaxis()->SetTitle("L1 Efficiency");
  histDummy->GetXaxis()->SetTitleSize(0.06); // default is 0.03                                                                    
  /* Set y-axis limits */
  histDummy->GetYaxis()->SetRangeUser(0.0, 1.1);

  /* Customize legend */
  for (itGraph = graphs.begin(), itLabel = labels.begin();
       itGraph != graphs.end();
       itGraph++, itLabel++)
    {
      leg->AddEntry(*itGraph, *itLabel,  "P");
    }

  gStyle->SetLegendFont(25);
  leg->Draw();

  Tcan->cd();
  Tcan->SaveAs(outputDir+outputName);
}
             

/*******************************************************************/

#endif