/*******************************************************************/
/* calculateEfficiency.cpp                                         */
/* Helper function                                                 */
/* Based on comparisonplots.C by Isobel Ojalvo                     */
/* Author: Stephanie Kwan                                          */
/*******************************************************************/

#include "../baseCodeForPlots/CMS_lumi.h"
#include "../baseCodeForPlots/tdrstyle.C"

#include "resolutionHist.cpp"

#include "TAxis.h"
#include "TChain.h"
#include "TF1.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TH1F.h"

#include "THStack.h"
#include "TLegend.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TTree.h"

#include <cmath>
#include <iostream>
#include <math.h>
#include <string>
#include <TMath.h>
#include <vector>

/* TMVA */
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"


#ifndef CALCULATE_RESOLUTION_CPP_INCL
#define CALCULATE_RESOLUTION_CPP_INCL

/*******************************************************************/

/* Helper function declarations */

//void setMaxErrorTo1(TGraphAsymmErrors *graph);

/*******************************************************************/

/* Calculates and returns the efficiency and statistical uncertainty
   for variable, using the n-tuple specified by rootFileDirectory
   and treePath, and applying the L1 cuts l1Cut and reco-level cut
   recoCut. Returns a TGraphAsymmErrors with x-axis range [low, high]).
   If variableBin = true, uses the bins specified in the function.
*/
   
TH1F* calculateResolution(TString variable,
                   TString treePath, TString rootFileDirectory,
                   TString l1Cut,
                   TString recoCut,
                   double low,
                   double high,
                   bool variableBin = false)
{
  /* Load file */
  TFile *file = new TFile(rootFileDirectory);
  if (!file->IsOpen() || file==0 )
    {
      std::cout<<"ERROR: FILE "<< rootFileDirectory <<" NOT FOUND; EXITING"<<std::endl;
      return NULL;
    }

  TTree* tree = (TTree*)file->Get(treePath);
  if (tree == 0)
    {
      std::cout<<"ERROR: Tree is "<< tree<<" NOT FOUND; EXITING"<<std::endl;
      return NULL;
    }


  /* Numerator and denominator histograms. */
  //TH1F* Num;
  int bins = 100;
  TH1F* Num = new TH1F("Num", "Num", bins, low, high);
  // Float_t xbins[10] = {0, 5, 10, 20, 25, 30, 50, 70, 100, 200};
  // Float_t xbins[11] = {20, 25, 30, 35, 40, 45, 50, 60, 70, 90, 110};
  //Float_t xbins[11] = {20, 25, 30, 35, 40, 45, 50, 60, 70, 90, 110};
  //int nVarBins = (int) sizeof(xbins)/sizeof(xbins[0]) - 1;

  //if(variableBin)
  //  {
  //    Denom = new TH1F("Denom", "Denom", nVarBins, xbins);
  //    Num = new TH1F("Num", "Num", nVarBins, xbins);
  //  }
  //else
  //  {
  //    Denom = new TH1F("Denom", "Denom", bins, low, high);
  //    Num = new TH1F("Num", "Num", bins, low, high);
  //  }
  //Denom->Sumw2();
  Num->Sumw2();

  /* Fill the histograms. */
  //tree->Draw(variable+">>+Denom", recoCut);
  //tree->Draw(variable+">>+Num", l1Cut);
  tree->Draw(variable+">>+Num", recoCut);

  /*
  for (int i = 0; i < 10; i++)
    {
      printf("Num bin %d content is %f, with error %f\n", i,
        Num->GetBinContent(i),
             Num->GetBinError(i));
    }
  */

  //Num->Divide(Denom);

  //TGraphAsymmErrors* effAsym = new TGraphAsymmErrors(Num);

  //setMaxErrorTo1(effAsym);

  return Num;
}

/*******************************************************************/

/* Sets the maximum and minimum error of graph to be 1.0 and 0.0 
   respectively. */

//void setMaxErrorTo1(TGraphAsymmErrors *graph)
//{
//  for (int i = 0; i < graph->GetN(); i++)
//    {
//      Double_t errorY = graph->GetErrorY(i);
//      Double_t pointX, pointY;
//
//      if (graph->GetPoint(i, pointX, pointY) < 0)
//        printf("Error getting point\n");
//
//      Double_t errorUp = pointY + errorY;
//      Double_t errorLow = pointY - errorY;
//
//      if (errorUp > 1)
//        graph->SetPointEYhigh(i, 1 - pointY);
//      else if (errorLow < 0)
//        graph->SetPointEYlow(i, pointY);
//
//    }
//}



/*******************************************************************/

#endif
