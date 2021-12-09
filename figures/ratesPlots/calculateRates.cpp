/*******************************************************************/
/* calculateRates.cpp                                              */
/* Helper function                                                 */
/* Based on comparisonplots.C by Isobel Ojalvo                     */
/* Author: Stephanie Kwan                                          */
/*******************************************************************/

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1.h"
#include "TF1.h"
#include <math.h>
#include <TEfficiency.h>
#include <TMath.h>
#include <iostream>
#include <string>
#include <vector>

#include <cmath>

#include "../baseCodeForPlots/tdrstyle.C"
#include "../baseCodeForPlots/CMS_lumi.h"
#include "../baseCodeForPlots/helperFunctions.C"

#ifndef CALCULATE_RATES_INCL
#define CALCULATE_RATES_INCL

/*******************************************************************/

/* histPath: Path to histogram containing the objects passing the
                desired trigger working point. e.g. 
                "L1TauAnalyzerRates/l1TauIsoTight_pt" contains the
                pTs of the taus passing the Tight isolation MVA
                requirement.
   nEventsHistPath: Path to TH1F of nEvents (counter of # of events originally processed)
   rootFileDirectory: Directory to the ROOT file. 
   Returns a TH1F of the rate assuming the HL-LHC 40 MHz event rate. */
   
TH1F* calculateRates(TString histPath,
                     TString nEventsHistPath,
                     TString rootFileDirectory)
{
  /* Load file */
  TFile *file = new TFile(rootFileDirectory);
  if (!file->IsOpen() || file==0 )
    {
      std::cout<<"ERROR FILE "<< rootFileDirectory <<" NOT FOUND; EXITING"<<std::endl;
      return 0;
    }

  TH1F* hist = (TH1F*)file->Get(histPath);
  if (hist == 0)
    {
      std::cout << "ERROR: " << histPath << " not found; EXITING"<<std::endl;
      return 0;
    }


  int nBins = hist->GetSize();

  float xMin = hist->GetBinLowEdge(0);
  float xMax = (hist->GetBinLowEdge(nBins) + hist->GetBinWidth(nBins));

  //  printf("nBins = %i, xMin = %f, xMax = %f\n", nBins, xMin, xMax);

  TH1F* ratesHist = new TH1F("Rates", "Rates", nBins+2, xMin, xMax);
  ratesHist->Sumw2();

  /* Loop through bins in the Rates histogram, in reverse order. */
  int Sum = 0;

  for (int i = nBins; i > 0; i--)
    {
      std::cout << "Bin number " << i << " out of " << nBins << ": sum is " << Sum << std::endl;
      Sum += hist->GetBinContent(i);
      ratesHist->SetBinContent(i, Sum);
    }

  /* Calculate (# of all events passing the BDT) / (# all events) */
  double nPass = hist->GetEntries();
  double nEvents = getEvents(rootFileDirectory, nEventsHistPath);

  /* Convert each bin to a fraction of total events. */
  float firstBin = ratesHist->GetBinContent(1);
  ratesHist->Scale((double) 1.00 / firstBin);
  
  /* kHz */
  ratesHist->Scale(40.0 * 1000000.0 / 1000.0);

  return ratesHist;
}

#endif