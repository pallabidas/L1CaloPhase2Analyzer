/*******************************************************************/
/* makeEfficienciesPlot.cpp                                        */
/* ROOT macro                                                      */
/* Usage: root -l -b -q makeEfficienciesPlot.cpp                   */
/*******************************************************************/

#include "efficiencyHist.cpp"
#include "calculateEfficiency.cpp"

#include <string>
/*********************************************************************/

/* Plots L1 RCT and GCT EGamma efficiency as a function of gen-level variables. */

void makeEfficienciesPlot(void)
{
  gROOT->ProcessLine(".L calculateEfficiency.cpp");

  /* Load the TTree. */
  TString treePath = "l1NtupleProducer/efficiencyTree";
  //  TString rootFileDirectory = "../test/analyzer.root";
  TString rootFileDirectory = "/afs/cern.ch/work/p/pdas/emulator_phase2/try/CMSSW_11_1_7/src/L1Trigger/L1CaloPhase2Analyzer/test/analyzer.root";
  TString outputDirectory = "/afs/cern.ch/work/p/pdas/www/emulator_phase2/";

  float xMin, xMax;
  TString genCut, l1Cut;
  bool useVariableBinning;

  std::vector<TGraphAsymmErrors*> vGraphs;
  std::vector<TString> vLabels;
  std::vector<int> vColors;

  /*******************************************************/
  /* efficiency as a function of genPt                  */
  /*******************************************************/

  // xMin = 0;
  // xMax = 100;
  // genCut  = "(abs(genEta) < 1.4841)";
  // l1Cut   = "(abs(genEta) < 1.4841) && (rct_cPt > 25)";
  // useVariableBinning = false;

  // TGraphAsymmErrors* loose = calculateEfficiency("genPt", treePath, rootFileDirectory,
  //               l1Cut + "&& rct_cPt > 30",
  //               genCut, xMin, xMax, useVariableBinning);
  // vGraphs.push_back(loose);
  // vLabels.push_back("L1 p_{T} > 30");
  // vColors.push_back(kGreen+2);

  // TGraphAsymmErrors* medium = calculateEfficiency("genPt", treePath, rootFileDirectory,
  //                l1Cut + "&& rct_cPt > 35",
  //                genCut, xMin, xMax, useVariableBinning);
  // vGraphs.push_back(medium);
  // vLabels.push_back("L1 p_{T} > 35");
  // vColors.push_back(kBlue);
  
  // TGraphAsymmErrors* tight = calculateEfficiency("genPt", treePath, rootFileDirectory,
  //               l1Cut + "&& rct_cPt > 40",
  //               genCut, xMin, xMax, useVariableBinning);
  // vGraphs.push_back(tight);
  // vLabels.push_back("L1 p_{T} > 40");
  // vColors.push_back(kBlack);

  // TGraphAsymmErrors* all = calculateEfficiency("genPt", treePath, rootFileDirectory,
  //             l1Cut,
  //             genCut, xMin, xMax, useVariableBinning);
  // vGraphs.push_back(all);
  // vLabels.push_back("No additional cut");
  // vColors.push_back(kRed);

  // plotNEfficiencies(vGraphs, vLabels, vColors,
  //       "Gen Electron p_{T} [GeV]",
  //       "Phase 2 RCT",                                                                
  //       "efficiency_genPt_barrel_RCT",        
  //       outputDirectory);    
  
  /*******************************************************/
  /* efficiency as a function of recoEta                 */
  /*******************************************************/
  // vGraphs.clear();  vLabels.clear();  vColors.clear();

  // xMin = -1.5;
  // xMax = 1.5;

  // genCut  = "(abs(genEta) < 1.4841)";
  // l1Cut   = "(abs(genEta) < 1.4841) && (rct_cPt > 25)";
  // useVariableBinning = false;

  // TGraphAsymmErrors* loose2 = calculateEfficiency("genEta", treePath, rootFileDirectory,
  //               l1Cut + "&& rct_cPt > 30",
  //               genCut, xMin, xMax, useVariableBinning);
  // vGraphs.push_back(loose2);
  // vLabels.push_back("L1 p_{T} > 30");
  // vColors.push_back(kGreen+2);

  // TGraphAsymmErrors* medium2 = calculateEfficiency("genEta", treePath, rootFileDirectory,
  //                l1Cut + "&& rct_cPt > 35",
  //                genCut, xMin, xMax, useVariableBinning);
  // vGraphs.push_back(medium2);
  // vLabels.push_back("L1 p_{T} > 35");
  // vColors.push_back(kBlue);
  
  // TGraphAsymmErrors* tight2 = calculateEfficiency("genEta", treePath, rootFileDirectory,
  //               l1Cut + "&& rct_cPt > 40",
  //               genCut, xMin, xMax, useVariableBinning);
  // vGraphs.push_back(tight2);
  // vLabels.push_back("L1 p_{T} > 40");
  // vColors.push_back(kBlack);

  // TGraphAsymmErrors* all2 = calculateEfficiency("genEta", treePath, rootFileDirectory,
  //             l1Cut,
  //             genCut, xMin, xMax, useVariableBinning);
  // vGraphs.push_back(all2);
  // vLabels.push_back("No additional cut");
  // vColors.push_back(kRed);

  // plotNEfficiencies(vGraphs, vLabels, vColors,
  //       "Gen Electron #eta",
  //       "Phase 2 RCT",                                                                
  //       "efficiency_genEta_barrel_RCT",        
  //       outputDirectory);    



  /*******************************************************/
  /* efficiency as a function of genPt: GCT              */
  /*******************************************************/

  vGraphs.clear();  vLabels.clear();  vColors.clear();
  xMin = 0;
  xMax = 100;
  genCut  = "(abs(genEta) < 1.4841)";
  l1Cut   = "(abs(genEta) < 1.4841) && (gct_cPt > 25)";
  useVariableBinning = false;

  TGraphAsymmErrors* loose3 = calculateEfficiency("genPt", treePath, rootFileDirectory,
                l1Cut + "&& gct_cPt > 30",
                genCut, xMin, xMax, useVariableBinning);
  vGraphs.push_back(loose3);
  vLabels.push_back("L1 p_{T} > 30");
  vColors.push_back(kGreen + 2);

  TGraphAsymmErrors* medium3 = calculateEfficiency("genPt", treePath, rootFileDirectory,
                 l1Cut + "&& gct_cPt > 35",
                 genCut, xMin, xMax, useVariableBinning);
  vGraphs.push_back(medium3);
  vLabels.push_back("L1 p_{T} > 35");
  vColors.push_back(kBlue);
  
  TGraphAsymmErrors* tight3 = calculateEfficiency("genPt", treePath, rootFileDirectory,
                l1Cut + "&& gct_cPt > 40",
                genCut, xMin, xMax, useVariableBinning);
  vGraphs.push_back(tight3);
  vLabels.push_back("L1 p_{T} > 40");
  vColors.push_back(kBlack);

  TGraphAsymmErrors* all3 = calculateEfficiency("genPt", treePath, rootFileDirectory,
              l1Cut,
              genCut, xMin, xMax, useVariableBinning);
  vGraphs.push_back(all3);
  vLabels.push_back("No additional cut");
  vColors.push_back(kRed);

  plotNEfficiencies(vGraphs, vLabels, vColors,
        "Gen Electron p_{T} [GeV]",
        "Phase 2 GCT",                                                                
        "efficiency_genPt_barrel_GCT",        
        outputDirectory);    
  

  /*******************************************************/
  /* efficiency as a function of genEta: GCT             */
  /*******************************************************/

  vGraphs.clear();  vLabels.clear();  vColors.clear();
  xMin = -1.5;
  xMax = 1.5;

  genCut  = "(abs(genEta) < 1.4841)";
  l1Cut   = "(abs(genEta) < 1.4841) && (gct_cPt > 25)";
  useVariableBinning = false;


  TGraphAsymmErrors* loose4 = calculateEfficiency("genEta", treePath, rootFileDirectory,
                l1Cut + "&& gct_cPt > 30",
                genCut, xMin, xMax, useVariableBinning);
  vGraphs.push_back(loose4);
  vLabels.push_back("L1 p_{T} > 30");
  vColors.push_back(kGreen+2);

  TGraphAsymmErrors* medium4 = calculateEfficiency("genEta", treePath, rootFileDirectory,
                 l1Cut + "&& gct_cPt > 35",
                 genCut, xMin, xMax, useVariableBinning);
  vGraphs.push_back(medium4);
  vLabels.push_back("L1 p_{T} > 35");
  vColors.push_back(kBlue);
  
  TGraphAsymmErrors* tight4 = calculateEfficiency("genEta", treePath, rootFileDirectory,
                l1Cut + "&& gct_cPt > 40",
                genCut, xMin, xMax, useVariableBinning);
  vGraphs.push_back(tight4);
  vLabels.push_back("L1 p_{T} > 40");
  vColors.push_back(kBlack);

  TGraphAsymmErrors* all4 = calculateEfficiency("genEta", treePath, rootFileDirectory,
              l1Cut,
              genCut, xMin, xMax, useVariableBinning);
  vGraphs.push_back(all4);
  vLabels.push_back("No additional cut");
  vColors.push_back(kRed);

  plotNEfficiencies(vGraphs, vLabels, vColors,
        "Gen Electron #eta",
        "Phase 2 GCT",                                                                
        "efficiency_genEta_barrel_GCT",        
        outputDirectory);    
  

  /*******************************************************/
  /* efficiency as a function of genPhi: GCT             */
  /*******************************************************/

  vGraphs.clear();  vLabels.clear();  vColors.clear();
  xMin = -3.142;
  xMax = 3.142;

  genCut  = "(abs(genEta) < 1.4841)";
  l1Cut   = "(abs(genEta) < 1.4841) && (gct_cPt > 25)";
  useVariableBinning = false;


  TGraphAsymmErrors* loose5 = calculateEfficiency("genPhi", treePath, rootFileDirectory,
                l1Cut + "&& gct_cPt > 30",
                genCut, xMin, xMax, useVariableBinning);
  vGraphs.push_back(loose5);
  vLabels.push_back("L1 p_{T} > 30");
  vColors.push_back(kGreen+2);

  TGraphAsymmErrors* medium5 = calculateEfficiency("genPhi", treePath, rootFileDirectory,
                 l1Cut + "&& gct_cPt > 35",
                 genCut, xMin, xMax, useVariableBinning);
  vGraphs.push_back(medium5);
  vLabels.push_back("L1 p_{T} > 35");
  vColors.push_back(kBlue);

  TGraphAsymmErrors* tight5 = calculateEfficiency("genPhi", treePath, rootFileDirectory,
                l1Cut + "&& gct_cPt > 40",
                genCut, xMin, xMax, useVariableBinning);
  vGraphs.push_back(tight5);
  vLabels.push_back("L1 p_{T} > 40");
  vColors.push_back(kBlack);

  TGraphAsymmErrors* all5 = calculateEfficiency("genPhi", treePath, rootFileDirectory,
              l1Cut,
              genCut, xMin, xMax, useVariableBinning);
  vGraphs.push_back(all5);
  vLabels.push_back("No additional cut");
  vColors.push_back(kRed);

  plotNEfficiencies(vGraphs, vLabels, vColors,
        "Gen Electron #phi",
        "Phase 2 GCT",
        "efficiency_genPhi_barrel_GCT",
        outputDirectory);

  
}
/*********************************************************************/
