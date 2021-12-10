/*******************************************************************/
/* makeRatesPlots.cpp                                              */
/* ROOT macro                                                      */
/* Usage: root -l -b -q makeRatesPlots.cpp                         */
/*******************************************************************/

#include "calculateRates.cpp"
#include "plotNRates.cpp"
#include "../baseCodeForPlots/CMS_lumi.h"
#include "../baseCodeForPlots/tdrstyle.C"

void plotHists(TH1F* h1, TString h1Label, 
          TH1F* h2, TString h2Label,
          TH1F* h3, TString h3Label,
          TString filename,
               TString outputDir);

void plotFiveRates(TH1F* h1, TString h1Label, int c1,
         TH1F* h2, TString h2Label, int c2,
         TH1F* h3, TString h3Label, int c3,
         TH1F* h4, TString h4Label, int c4,
         TH1F* h5, TString h5Label, int c5,
         float xMin, float xMax,
         float yMin, float yMax,
         TString legendTitle,
         TString filename,
         TString outputDir);

/*********************************************************************/

// Create rates plots.

void makeRatesPlots(void)
{
  gROOT->ProcessLine(".L calculateRates.cpp");

  /* Specify paths/directories for input files. */
  TString rootFileDirectory = "/Users/stephaniekwan/Dropbox/Princeton_G4/Phase2RCT/analyzer/analyzer-rates.root";
  TString outputDirectory = "/Users/stephaniekwan/Dropbox/Princeton_G4/Phase2RCT/analyzer/ratesPlots/";

  float xMin, xMax;
  float yMin, yMax;
  std::vector<TH1F*> vHists;
  std::vector<TString> vLabels;
  std::vector<int> vColors;
  bool useLogy;


  /*******************************************************/
  /* Rates as a function of l1Pt                         */
  /*******************************************************/

  xMin = 0;
  xMax = 60.0;
  yMin = 10.0;
  yMax = 50000;
  useLogy = true;

  /* EG rates in barrel */
  TString folder= "l1NtupleProducer/";
  TString evHist = "l1NtupleProducer/nEvents";
  TH1F* egVLoose = calculateRates(folder + "l1eg_pt", evHist, rootFileDirectory);
  // TH1F* egLoose  = calculateRates(folder + "l1egLoose_pt",  evHist, rootFileDirectory);
  // TH1F* egMedium = calculateRates(folder + "l1egMedium_pt", evHist, rootFileDirectory);
  // TH1F* egTight  = calculateRates(folder + "l1egTight_pt",  evHist, rootFileDirectory);

  vHists.push_back(egVLoose); vLabels.push_back("No additional cuts"); vColors.push_back(kPink);
  // vHists.push_back(egLoose);  vLabels.push_back("L1 p_{T} > 30"); vColors.push_back(kTeal-8);
  // vHists.push_back(egMedium); vLabels.push_back("L1 p_{T} > 35"); vColors.push_back(kAzure-9);
  // vHists.push_back(egTight);  vLabels.push_back("L1 p_{T} > 40"); vColors.push_back(kAzure+2);
  
  plotNRates(vHists, vLabels, vColors,
             xMin, xMax, yMin, yMax,
             "L1 EGamma p_{T}",
             "rates_barrel_GCT",
             outputDirectory,
             useLogy);

  vHists.clear();  vLabels.clear();  vColors.clear();

}

// /*********************************************************************/

// /* Plot five histograms with rates, with the labels h1Label and line
//    color int c1 etc. If there's a way to do this without hard-coding
//    five histograms in the argument, I'd love to see it... */

// void plotFiveRates(TH1F* h1, TString h1Label, int c1,
//          TH1F* h2, TString h2Label, int c2,
//          TH1F* h3, TString h3Label, int c3,
//          TH1F* h4, TString h4Label, int c4,
//          TH1F* h5, TString h5Label, int c5,
//          float xMin, float xMax,
//          float yMin, float yMax,
//          TString legendTitle,
//          TString filename,
//          TString outputDir)
// {
//   /*******************************************************/
//   /* plotting                                            */
//   /*******************************************************/
//   setTDRStyle();
//   TCanvas* Tcan = new TCanvas("Tcan","", 100, 20, 1000, 1000);
//   TLegend* leg = new TLegend(0.50, 0.6, 0.90, 0.90);
//   Tcan->SetGrid();

//   Tcan->cd();     /* Set current canvas */
//   Tcan->SetFillColor(0);
//   applyLegStyle(leg);

//   h1->SetLineWidth(3);
//   h1->SetLineColor(c1);

//   h2->SetLineWidth(3);
//   h2->SetLineColor(c2);

//   h3->SetLineWidth(3);
//   h3->SetLineColor(c3);

//   h4->SetLineWidth(3);
//   h4->SetLineColor(c4);

//   h5->SetLineWidth(3);
//   h5->SetLineColor(c5);

//   /* Set x-axis limits */
//   h1->Draw("");
//   h2->Draw("SAME");
//   h3->Draw("SAME");
//   h4->Draw("SAME");
//   h5->Draw("SAME");

//   h1->GetXaxis()->SetRangeUser(xMin, xMax);
//   h1->GetYaxis()->SetRangeUser(yMin, yMax);
//   h1->GetXaxis()->SetTitle("Level-1 #tau_{H} p_{T} [GeV]");
//   h1->GetYaxis()->SetTitle("Rate [kHz]");
//   h1->GetXaxis()->SetTitleSize(0.06); // default is 0.03                                                                                                        
//   h1->GetYaxis()->SetTitleSize(0.06);


//   Tcan->SetLogy();

//   /* Customize legend */
//   leg->SetHeader(legendTitle);
//   leg->AddEntry(h1, h1Label, "l");
//   leg->AddEntry(h2, h2Label, "l");
//   leg->AddEntry(h3, h3Label, "l");
//   leg->AddEntry(h4, h4Label, "l");
//   leg->AddEntry(h5, h5Label, "l");

//   gStyle->SetLegendFont(20);
//   leg->Draw();

//   Tcan->cd();
//   Tcan->SaveAs(outputDir+filename);

// }


// /*********************************************************************/

// /* Plots Hist. */

// void plotHists(TH1F* h1, TString h1Label,
//                TH1F* h2, TString h2Label,
//                TH1F* h3, TString h3Label,
//           TString filename,
//           TString outputDir)
// {
//   /*******************************************************/
//   /* plotting                                            */
//   /*******************************************************/
//   setTDRStyle();
//   TCanvas* Tcan = new TCanvas("Tcan","", 100, 20, 1000, 800);
//   TLegend* leg = new TLegend(0.60,0.75,0.85,0.9);
//   Tcan->SetGrid();

//   Tcan->cd();     /* Set current canvas */
//   Tcan->SetFillColor(0);
//   leg = new TLegend(0.55, 0.2, 0.90, 0.5);
//   applyLegStyle(leg);

//   h1->SetMarkerColor(kViolet-5);
//   h1->SetLineWidth(2);
//   h1->SetLineColor(kViolet-5);

//   h2->SetMarkerColor(kAzure+7);
//   h2->SetLineWidth(2);
//   h2->SetLineColor(kAzure+7);

//   h3->SetMarkerColor(kOrange+7);
//   h3->SetLineWidth(2);
//   h3->SetLineColor(kOrange+7);

//   /* Set x-axis limits */
//   h1->GetXaxis()->SetRangeUser(0.0, 90.0);
  
//   h1->Draw("");
//   h2->Draw("SAME");
//   h3->Draw("SAME");  

//   h1->GetXaxis()->SetTitle("L1 Tau p_{T} [GeV]");
//   h1->GetYaxis()->SetTitle("Rate [Hz]");
//   h1->GetXaxis()->SetTitleSize(0.06); // default is 0.03  
//   h1->GetYaxis()->SetTitleSize(0.06);

  
//   Tcan->SetLogy();

//   /* Customize legend */
//   leg->SetHeader("Phase 2 L1 Taus");
//   leg->AddEntry(h1, h1Label, "l");
//   leg->AddEntry(h2, h2Label, "l");
//   leg->AddEntry(h3, h3Label, "l");

//   gStyle->SetLegendFont(20);
//   leg->Draw();

//   Tcan->cd();
//   Tcan->SaveAs(outputDir+filename);

// }

/*********************************************************************/