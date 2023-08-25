/*******************************************************************/
/* makeEfficienciesPlot.cpp                                        */
/* ROOT macro                                                      */
/* Usage: root -l -b -q makeEfficienciesPlot_pion.cpp              */
/*******************************************************************/

#include "efficiencyHist.cpp"
#include "calculateEfficiency.cpp"
#include "resolutionHist.cpp"
#include "calculateResolution.cpp"

#include <string>
/*********************************************************************/


void makeEfficienciesPlotJet(void)
{
  gROOT->ProcessLine(".L calculateEfficiency.cpp");

  /* Load the TTree. */
  TString treePath = "l1NtupleProducer/jetEfficiencyTree";
  TString rootFileDirectory = "/afs/cern.ch/work/p/pdas/emulator_phase2/calojet/CMSSW_13_2_0/src/L1Trigger/L1CaloPhase2Analyzer/test/analyzer_VBFH.root";
  TString outputDirectory = "/afs/cern.ch/work/p/pdas/www/emulator_phase2/13_2_0/";

  float xMin, xMax;
  TString genJetCut, l1Cut;
  bool useVariableBinning;

  std::vector<TGraphAsymmErrors*> vGraphs;
  std::vector<TString> vLabels;
  std::vector<int> vColors;
  std::vector<TH1F*> vReso;
  std::vector<TString> vResoLabels;
  std::vector<int> vResoColors;

  /*******************************************************/
  /* resolution as a function of genJetPt                */
  /*******************************************************/

  vReso.clear(); vResoLabels.clear(); vResoColors.clear();

  TH1F* reso1 = calculateResolution("(gctJet_Pt - genJetPt)/genJetPt", treePath, rootFileDirectory,
               "(abs(genJetEta) < 1.4841) && genJetPt > 30 && gctJet_deltaR < 0.4", -1, 5, useVariableBinning);
  vReso.push_back(reso1);
  vResoLabels.push_back("barrel");
  vResoColors.push_back(kRed);

  TH1F* reso2 = calculateResolution("(gctJet_Pt - genJetPt)/genJetPt", treePath, rootFileDirectory,
               "(abs(genJetEta) > 1.4841) && (abs(genJetEta) < 3.0) && genJetPt > 30 && gctJet_deltaR < 0.4", -1, 5, useVariableBinning);
  vReso.push_back(reso2);
  vResoLabels.push_back("endcap");
  vResoColors.push_back(kAzure+1);

  TH1F* reso3 = calculateResolution("(gctJet_Pt - genJetPt)/genJetPt", treePath, rootFileDirectory,
               "(abs(genJetEta) > 3.0) && (abs(genJetEta) < 5.0) && genJetPt > 30 && gctJet_deltaR < 0.4", -1, 5, useVariableBinning);
  vReso.push_back(reso3);
  vResoLabels.push_back("forward");
  vResoColors.push_back(kGreen+1);

  plotNResolutions(vReso, vResoLabels, vResoColors,
        "Resolution vs Gen p_{T}",
        "resolution_genJetPt",
        outputDirectory);

  /*******************************************************/
  /* resolution as a function of genJetEta               */
  /*******************************************************/

  vReso.clear(); vResoLabels.clear(); vResoColors.clear();

  TH1F* reso4 = calculateResolution("(gctJet_Eta - genJetEta)/genJetEta", treePath, rootFileDirectory,
               "(abs(genJetEta) < 1.4841) && genJetPt > 30 && gctJet_deltaR < 0.4", -1, 1, useVariableBinning);
  vReso.push_back(reso4);
  vResoLabels.push_back("barrel");
  vResoColors.push_back(kRed);

  TH1F* reso5 = calculateResolution("(gctJet_Eta - genJetEta)/genJetEta", treePath, rootFileDirectory,
               "(abs(genJetEta) > 1.4841) && (abs(genJetEta) < 3.0) && genJetPt > 30 && gctJet_deltaR < 0.4", -1, 1, useVariableBinning);
  vReso.push_back(reso5);
  vResoLabels.push_back("endcap");
  vResoColors.push_back(kAzure+1);

  TH1F* reso6 = calculateResolution("(gctJet_Eta - genJetEta)/genJetEta", treePath, rootFileDirectory,
               "(abs(genJetEta) > 3.0) && (abs(genJetEta) < 5.0) && genJetPt > 30 && gctJet_deltaR < 0.4", -1, 1, useVariableBinning);
  vReso.push_back(reso6);
  vResoLabels.push_back("forward");
  vResoColors.push_back(kGreen+1);

  plotNResolutions(vReso, vResoLabels, vResoColors,
        "Resolution vs Gen #eta",
        "resolution_genJetEta",
        outputDirectory);

  /*******************************************************/
  /* resolution as a function of genJetPhi               */
  /*******************************************************/

  vReso.clear(); vResoLabels.clear(); vResoColors.clear();

  TH1F* reso7 = calculateResolution("(gctJet_Phi - genJetPhi)/genJetPhi", treePath, rootFileDirectory,
               "(abs(genJetEta) < 1.4841) && genJetPt > 30 && gctJet_deltaR < 0.4", -1, 1, useVariableBinning);
  vReso.push_back(reso7);
  vResoLabels.push_back("barrel");
  vResoColors.push_back(kRed);

  TH1F* reso8 = calculateResolution("(gctJet_Phi - genJetPhi)/genJetPhi", treePath, rootFileDirectory,
               "(abs(genJetEta) > 1.4841) && (abs(genJetEta) < 3.0) && genJetPt > 30 && gctJet_deltaR < 0.4", -1, 1, useVariableBinning);
  vReso.push_back(reso8);
  vResoLabels.push_back("endcap");
  vResoColors.push_back(kAzure+1);

  TH1F* reso9 = calculateResolution("(gctJet_Phi - genJetPhi)/genJetPhi", treePath, rootFileDirectory,
               "(abs(genJetEta) > 3.0) && (abs(genJetEta) < 5.0) && genJetPt > 30 && gctJet_deltaR < 0.4", -1, 1, useVariableBinning);
  vReso.push_back(reso9);
  vResoLabels.push_back("forward");
  vResoColors.push_back(kGreen+1);

  plotNResolutions(vReso, vResoLabels, vResoColors,
        "Resolution vs Gen #phi",
        "resolution_genJetPhi",
        outputDirectory);

  /*******************************************************/
  /* efficiency as a function of genJetPt                */
  /*******************************************************/

  vGraphs.clear();  vLabels.clear();  vColors.clear();
  xMin = 0;
  xMax = 250;
  //genJetCut  = "(abs(genJetEta) < 5.0)";
  //l1Cut   = "(abs(genJetEta) < 1.4841) && gctJet_deltaR < 0.4";
  useVariableBinning = false;

  TGraphAsymmErrors* loose1 = calculateEfficiency("genJetPt", treePath, rootFileDirectory,
               "gctJet_deltaR < 0.4 && (abs(genJetEta) < 1.4841) && genJetPt > 0",
               "(abs(genJetEta) < 1.4841) && genJetPt > 0", xMin, xMax, useVariableBinning);
  vGraphs.push_back(loose1);
  vLabels.push_back("p_{T} > 0");
  vColors.push_back(kGreen + 2);

  TGraphAsymmErrors* medium1 = calculateEfficiency("genJetPt", treePath, rootFileDirectory,
               "gctJet_deltaR < 0.4 && (abs(genJetEta) < 1.4841) && genJetPt > 10",
               "(abs(genJetEta) < 1.4841) && genJetPt > 10", xMin, xMax, useVariableBinning);
  vGraphs.push_back(medium1);
  vLabels.push_back("p_{T} > 10");
  vColors.push_back(kBlue);
  
  TGraphAsymmErrors* tight1 = calculateEfficiency("genJetPt", treePath, rootFileDirectory,
               "gctJet_deltaR < 0.4 && (abs(genJetEta) < 1.4841) && genJetPt > 30",
               "(abs(genJetEta) < 1.4841) && genJetPt > 30", xMin, xMax, useVariableBinning);
  vGraphs.push_back(tight1);
  vLabels.push_back("p_{T} > 30");
  vColors.push_back(kBlack);

  TGraphAsymmErrors* all1 = calculateEfficiency("genJetPt", treePath, rootFileDirectory,
               "gctJet_deltaR < 0.4 && (abs(genJetEta) < 1.4841) && genJetPt > 100",
               "(abs(genJetEta) < 1.4841) && genJetPt > 100", xMin, xMax, useVariableBinning);
  vGraphs.push_back(all1);
  vLabels.push_back("p_{T} > 100");
  vColors.push_back(kRed);

  plotNEfficiencies(vGraphs, vLabels, vColors,
        "Gen p_{T} [GeV]",
        "Phase 2 Calo",                                                                
        "efficiency_genJetPt_barrel",
        outputDirectory);    

  vGraphs.clear();  vLabels.clear();  vColors.clear();
  vGraphs.push_back(all1);
  vLabels.push_back("barrel");
  vColors.push_back(kRed);

  TGraphAsymmErrors* all2 = calculateEfficiency("genJetPt", treePath, rootFileDirectory,
               "gctJet_deltaR < 0.4 && (abs(genJetEta) > 1.4841) && (abs(genJetEta) < 3.0) && genJetPt > 100",
               "(abs(genJetEta) > 1.4841) && (abs(genJetEta) < 3.0) && genJetPt > 100", xMin, xMax, useVariableBinning);
  vGraphs.push_back(all2);
  vLabels.push_back("endcap");
  vColors.push_back(kAzure+1);

  TGraphAsymmErrors* all3 = calculateEfficiency("genJetPt", treePath, rootFileDirectory,
               "gctJet_deltaR < 0.4 && (abs(genJetEta) > 3.0) && (abs(genJetEta) < 5.0) && genJetPt > 100",
               "(abs(genJetEta) > 3.0) && (abs(genJetEta) < 5.0) && genJetPt > 100", xMin, xMax, useVariableBinning);
  vGraphs.push_back(all3);
  vLabels.push_back("forward");
  vColors.push_back(kGreen+1);

  plotNEfficiencies(vGraphs, vLabels, vColors,
        "Gen p_{T} [GeV]",
        "Phase 2 Calo",
        "efficiency_genJetPt",
        outputDirectory);
  
  /*******************************************************/
  /* efficiency as a function of genJetEta               */
  /*******************************************************/

  vGraphs.clear();  vLabels.clear();  vColors.clear();
  xMin = -1.4841;
  xMax = 1.4841;

  //genJetCut  = "(abs(genJetEta) < 1.4841)";
  //l1Cut   = "(abs(genJetEta) < 1.4841) && gctJet_deltaR < 0.4";
  useVariableBinning = false;

  TGraphAsymmErrors* loose2 = calculateEfficiency("genJetEta", treePath, rootFileDirectory,
               "gctJet_deltaR < 0.4 && (abs(genJetEta) < 1.4841) && genJetPt > 0",
               "(abs(genJetEta) < 1.4841) && genJetPt > 0", xMin, xMax, useVariableBinning);
  vGraphs.push_back(loose2);
  vLabels.push_back("p_{T} > 0");
  vColors.push_back(kGreen+2);

  TGraphAsymmErrors* medium2 = calculateEfficiency("genJetEta", treePath, rootFileDirectory,
               "gctJet_deltaR < 0.4 && (abs(genJetEta) < 1.4841) && genJetPt > 10",
               "(abs(genJetEta) < 1.4841) && genJetPt > 10", xMin, xMax, useVariableBinning);
  vGraphs.push_back(medium2);
  vLabels.push_back("p_{T} > 10");
  vColors.push_back(kBlue);
  
  TGraphAsymmErrors* tight2 = calculateEfficiency("genJetEta", treePath, rootFileDirectory,
               "gctJet_deltaR < 0.4 && (abs(genJetEta) < 1.4841) && genJetPt > 30",
               "(abs(genJetEta) < 1.4841) && genJetPt > 30", xMin, xMax, useVariableBinning);
  vGraphs.push_back(tight2);
  vLabels.push_back("p_{T} > 30");
  vColors.push_back(kBlack);

  TGraphAsymmErrors* all4 = calculateEfficiency("genJetEta", treePath, rootFileDirectory,
               "gctJet_deltaR < 0.4 && (abs(genJetEta) < 1.4841) && genJetPt > 100",
               "(abs(genJetEta) < 1.4841) && genJetPt > 100", xMin, xMax, useVariableBinning);
  vGraphs.push_back(all4);
  vLabels.push_back("p_{T} > 100");
  vColors.push_back(kRed);

  plotNEfficiencies(vGraphs, vLabels, vColors,
        "Gen #eta",
        "Phase 2 Calo",                                                                
        "efficiency_genJetEta_barrel",
        outputDirectory);    

  vGraphs.clear();  vLabels.clear();  vColors.clear();
  xMin = -5.;
  xMax = 5.;

  TGraphAsymmErrors* all5 = calculateEfficiency("genJetEta", treePath, rootFileDirectory,
               "gctJet_deltaR < 0.4 && (abs(genJetEta) < 1.4841) && genJetPt > 30",
               "(abs(genJetEta) < 1.4841) && genJetPt > 30", xMin, xMax, useVariableBinning);
  vGraphs.push_back(all5);
  vLabels.push_back("barrel");
  vColors.push_back(kRed);

  TGraphAsymmErrors* all6 = calculateEfficiency("genJetEta", treePath, rootFileDirectory,
               "gctJet_deltaR < 0.4 && (abs(genJetEta) > 1.4841) && (abs(genJetEta) < 3.0) && genJetPt > 30",
               "(abs(genJetEta) > 1.4841) && (abs(genJetEta) < 3.0) && genJetPt > 30", xMin, xMax, useVariableBinning);
  vGraphs.push_back(all6);
  vLabels.push_back("endcap");
  vColors.push_back(kAzure+1);

  TGraphAsymmErrors* all7 = calculateEfficiency("genJetEta", treePath, rootFileDirectory,
               "gctJet_deltaR < 0.4 && (abs(genJetEta) > 3.0) && (abs(genJetEta) < 5.0) && genJetPt > 30",
               "(abs(genJetEta) > 3.0) && (abs(genJetEta) < 5.0) && genJetPt > 30", xMin, xMax, useVariableBinning);
  vGraphs.push_back(all7);
  vLabels.push_back("forward");
  vColors.push_back(kGreen+1);

  plotNEfficiencies(vGraphs, vLabels, vColors,
        "Gen #eta",
        "Phase 2 Calo",
        "efficiency_genJetEta",
        outputDirectory);

  
  /*******************************************************/
  /* efficiency as a function of genJetPhi               */
  /*******************************************************/

  vGraphs.clear();  vLabels.clear();  vColors.clear();
  xMin = -3.142;
  xMax = 3.142;

  //genJetCut  = "(abs(genJetEta) < 1.4841)";
  //l1Cut   = "(abs(genJetEta) < 1.4841) && gctJet_deltaR < 0.4";
  useVariableBinning = false;


  TGraphAsymmErrors* loose3 = calculateEfficiency("genJetPhi", treePath, rootFileDirectory,
               "gctJet_deltaR < 0.4 && (abs(genJetEta) > 3.0) && (abs(genJetEta) < 5.0) && genJetPt > 0",
               "(abs(genJetEta) > 3.0) && (abs(genJetEta) < 5.0) && genJetPt > 0", xMin, xMax, useVariableBinning);
  vGraphs.push_back(loose3);
  vLabels.push_back("p_{T} > 0");
  vColors.push_back(kGreen+2);

  TGraphAsymmErrors* medium3 = calculateEfficiency("genJetPhi", treePath, rootFileDirectory,
               "gctJet_deltaR < 0.4 && (abs(genJetEta) > 3.0) && (abs(genJetEta) < 5.0) && genJetPt > 10",
               "(abs(genJetEta) > 3.0) && (abs(genJetEta) < 5.0) && genJetPt > 10", xMin, xMax, useVariableBinning);
  vGraphs.push_back(medium3);
  vLabels.push_back("p_{T} > 10");
  vColors.push_back(kBlue);

  TGraphAsymmErrors* tight3 = calculateEfficiency("genJetPhi", treePath, rootFileDirectory,
               "gctJet_deltaR < 0.4 && (abs(genJetEta) > 3.0) && (abs(genJetEta) < 5.0) && genJetPt > 30",
               "(abs(genJetEta) > 3.0) && (abs(genJetEta) < 5.0) && genJetPt > 30", xMin, xMax, useVariableBinning);
  vGraphs.push_back(tight3);
  vLabels.push_back("p_{T} > 30");
  vColors.push_back(kBlack);

  TGraphAsymmErrors* all8 = calculateEfficiency("genJetPhi", treePath, rootFileDirectory,
               "gctJet_deltaR < 0.4 && (abs(genJetEta) > 3.0) && (abs(genJetEta) < 5.0) && genJetPt > 100",
               "(abs(genJetEta) > 3.0) && (abs(genJetEta) < 5.0) && genJetPt > 100", xMin, xMax, useVariableBinning);
  vGraphs.push_back(all8);
  vLabels.push_back("p_{T} > 100");
  vColors.push_back(kRed);

  plotNEfficiencies(vGraphs, vLabels, vColors,
        "Gen #phi",
        "Phase 2 Calo",
        "efficiency_genJetPhi_barrel",
        outputDirectory);

  vGraphs.clear();  vLabels.clear();  vColors.clear();
  vGraphs.push_back(tight3);
  vLabels.push_back("barrel");
  vColors.push_back(kRed);

  TGraphAsymmErrors* all9 = calculateEfficiency("genJetPhi", treePath, rootFileDirectory,
               "gctJet_deltaR < 0.4 && (abs(genJetEta) > 1.4841) && (abs(genJetEta) < 3.0) && genJetPt > 30",
               "(abs(genJetEta) > 1.4841) && (abs(genJetEta) < 3.0) && genJetPt > 30", xMin, xMax, useVariableBinning);
  vGraphs.push_back(all9);
  vLabels.push_back("endcap");
  vColors.push_back(kAzure+1);

  TGraphAsymmErrors* all10 = calculateEfficiency("genJetPhi", treePath, rootFileDirectory,
               "gctJet_deltaR < 0.4 && (abs(genJetEta) > 3.0) && (abs(genJetEta) < 5.0) && genJetPt > 30",
               "(abs(genJetEta) > 3.0) && (abs(genJetEta) < 5.0) && genJetPt > 30", xMin, xMax, useVariableBinning);
  vGraphs.push_back(all10);
  vLabels.push_back("forward");
  vColors.push_back(kGreen+1);

  plotNEfficiencies(vGraphs, vLabels, vColors,
        "Gen #phi",
        "Phase 2 Calo",
        "efficiency_genJetPhi",
        outputDirectory);
  
}
/*********************************************************************/
