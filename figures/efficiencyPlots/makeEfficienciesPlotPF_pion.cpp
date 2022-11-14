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

/* Plots L1 RCT and PF cluster efficiency as a function of gen-level variables. */

void makeEfficienciesPlotPF_pion(void)
{
  gROOT->ProcessLine(".L calculateEfficiency.cpp");

  /* Load the TTree. */
  TString treePath = "l1NtupleProducer/efficiencyTree";
  TString rootFileDirectory = "/afs/cern.ch/work/p/pdas/emulator_phase2/CMSSW_12_3_0_pre4/src/L1Trigger/L1CaloPhase2Analyzer/test/analyzer_singlepion.root";
  TString outputDirectory = "/afs/cern.ch/work/p/pdas/www/emulator_phase2/12_3_0_pre4/singlepion/";

  float xMin, xMax;
  TString genCut, l1Cut;
  bool useVariableBinning;

  std::vector<TGraphAsymmErrors*> vGraphs;
  std::vector<TString> vLabels;
  std::vector<int> vColors;
  std::vector<TH1F*> vReso;
  std::vector<TString> vResoLabels;
  std::vector<int> vResoColors;

  /*******************************************************/
  /* resolution as a function of genPt  c                */
  /*******************************************************/

  vReso.clear(); vResoLabels.clear(); vResoColors.clear();

  TH1F* reso1 = calculateResolution("(pf_cPt - genPt)/genPt", treePath, rootFileDirectory,
                "(abs(genEta) < 1.4841) && genPt > 5 && pf_deltaR < 0.2", -1, 1, useVariableBinning);
  vReso.push_back(reso1);
  vResoLabels.push_back("PF cluster");
  vResoColors.push_back(kRed);

  TH1F* reso2 = calculateResolution("(gct_cPt - genPt)/genPt", treePath, rootFileDirectory,
                "(abs(genEta) < 1.4841) && genPt > 5 && gct_deltaR < 0.2", -1, 1, useVariableBinning);
  vReso.push_back(reso2);
  vResoLabels.push_back("EG cluster");
  vResoColors.push_back(kBlack);

  plotNResolutions(vReso, vResoLabels, vResoColors,
        "Resolution vs Gen p_{T}",
        "resolution_genPt_barrel",
        outputDirectory);

  /*******************************************************/
  /* resolution as a function of genEtac                 */
  /*******************************************************/

  vReso.clear(); vResoLabels.clear(); vResoColors.clear();

  TH1F* reso3 = calculateResolution("(pf_cEta - genEta)/genEta", treePath, rootFileDirectory,
                "(abs(genEta) < 1.4841) && genPt > 5 && pf_deltaR < 0.2", -1, 1, useVariableBinning);
  vReso.push_back(reso3);
  vResoLabels.push_back("PF cluster");
  vResoColors.push_back(kRed);

  TH1F* reso4 = calculateResolution("(gct_cEta - genEta)/genEta", treePath, rootFileDirectory,
                "(abs(genEta) < 1.4841) && genPt > 5 && gct_deltaR < 0.2", -1, 1, useVariableBinning);
  vReso.push_back(reso4);
  vResoLabels.push_back("EG cluster");
  vResoColors.push_back(kBlack);

  plotNResolutions(vReso, vResoLabels, vResoColors,
        "Resolution vs Gen #eta",
        "resolution_genEta_barrel",
        outputDirectory);

  /*******************************************************/
  /* resolution as a function of genEta                  */
  /*******************************************************/

  vReso.clear(); vResoLabels.clear(); vResoColors.clear();

  TH1F* reso5 = calculateResolution("(pf_cPhi - genPhi)/genPhi", treePath, rootFileDirectory,
                "(abs(genEta) < 1.4841) && genPt > 5 && pf_deltaR < 0.2", -1, 1, useVariableBinning);
  vReso.push_back(reso5);
  vResoLabels.push_back("PF cluster");
  vResoColors.push_back(kRed);

  TH1F* reso6 = calculateResolution("(gct_cPhi - genPhi)/genPhi", treePath, rootFileDirectory,
                "(abs(genEta) < 1.4841) && genPt > 5 && gct_deltaR < 0.2", -1, 1, useVariableBinning);
  vReso.push_back(reso6);
  vResoLabels.push_back("EG cluster");
  vResoColors.push_back(kBlack);

  plotNResolutions(vReso, vResoLabels, vResoColors,
        "Resolution vs Gen #phi",
        "resolution_genPhi_barrel",
        outputDirectory);

  /*******************************************************/
  /* efficiency as a function of genPt: PF               */
  /*******************************************************/

  vGraphs.clear();  vLabels.clear();  vColors.clear();
  xMin = 0;
  xMax = 200;
  genCut  = "(abs(genEta) < 1.4841)";
  l1Cut   = "(abs(genEta) < 1.4841) && pf_deltaR < 0.2";
  useVariableBinning = false;

  TGraphAsymmErrors* loose3 = calculateEfficiency("genPt", treePath, rootFileDirectory,
                l1Cut + "&& pf_cPt > 0",
                genCut, xMin, xMax, useVariableBinning);
  vGraphs.push_back(loose3);
  vLabels.push_back("p_{T} > 0");
  vColors.push_back(kGreen + 2);

  TGraphAsymmErrors* medium3 = calculateEfficiency("genPt", treePath, rootFileDirectory,
                 l1Cut + "&& pf_cPt > 2",
                 genCut, xMin, xMax, useVariableBinning);
  vGraphs.push_back(medium3);
  vLabels.push_back("p_{T} > 2");
  vColors.push_back(kBlue);
  
  TGraphAsymmErrors* tight3 = calculateEfficiency("genPt", treePath, rootFileDirectory,
                l1Cut + "&& pf_cPt > 5",
                genCut, xMin, xMax, useVariableBinning);
  vGraphs.push_back(tight3);
  vLabels.push_back("p_{T} > 5");
  vColors.push_back(kBlack);

  TGraphAsymmErrors* all3 = calculateEfficiency("genPt", treePath, rootFileDirectory,
              l1Cut + "&& pf_cPt > 10",
              genCut, xMin, xMax, useVariableBinning);
  vGraphs.push_back(all3);
  vLabels.push_back("p_{T} > 10");
  vColors.push_back(kRed);

  plotNEfficiencies(vGraphs, vLabels, vColors,
        "Gen p_{T} [GeV]",
        "Phase 2 GCT",                                                                
        "efficiency_genPt_barrel_PF_lowPt",        
        outputDirectory);    
  
  /*******************************************************/
  /* efficiency as a function of genEta: PF              */
  /*******************************************************/

  vGraphs.clear();  vLabels.clear();  vColors.clear();
  xMin = -1.5;
  xMax = 1.5;

  genCut  = "(abs(genEta) < 1.4841)";
  l1Cut   = "(abs(genEta) < 1.4841) && pf_deltaR < 0.2";
  useVariableBinning = false;


  TGraphAsymmErrors* loose4 = calculateEfficiency("genEta", treePath, rootFileDirectory,
                l1Cut + "&& pf_cPt > 0",
                genCut, xMin, xMax, useVariableBinning);
  vGraphs.push_back(loose4);
  vLabels.push_back("p_{T} > 0");
  vColors.push_back(kGreen+2);

  TGraphAsymmErrors* medium4 = calculateEfficiency("genEta", treePath, rootFileDirectory,
                 l1Cut + "&& pf_cPt > 2",
                 genCut, xMin, xMax, useVariableBinning);
  vGraphs.push_back(medium4);
  vLabels.push_back("p_{T} > 2");
  vColors.push_back(kBlue);
  
  TGraphAsymmErrors* tight4 = calculateEfficiency("genEta", treePath, rootFileDirectory,
                l1Cut + "&& pf_cPt > 5",
                genCut, xMin, xMax, useVariableBinning);
  vGraphs.push_back(tight4);
  vLabels.push_back("p_{T} > 5");
  vColors.push_back(kBlack);

  TGraphAsymmErrors* all4 = calculateEfficiency("genEta", treePath, rootFileDirectory,
              l1Cut + "&& pf_cPt > 10",
              genCut, xMin, xMax, useVariableBinning);
  vGraphs.push_back(all4);
  vLabels.push_back("p_{T} > 10");
  vColors.push_back(kRed);

  plotNEfficiencies(vGraphs, vLabels, vColors,
        "Gen #eta",
        "Phase 2 GCT",                                                                
        "efficiency_genEta_barrel_PF_lowPt",        
        outputDirectory);    
  
  /*******************************************************/
  /* efficiency as a function of genPhi: PF              */
  /*******************************************************/

  vGraphs.clear();  vLabels.clear();  vColors.clear();
  xMin = -3.142;
  xMax = 3.142;

  genCut  = "(abs(genEta) < 1.4841)";
  l1Cut   = "(abs(genEta) < 1.4841) && pf_deltaR < 0.2";
  useVariableBinning = false;


  TGraphAsymmErrors* loose5 = calculateEfficiency("genPhi", treePath, rootFileDirectory,
                l1Cut + "&& pf_cPt > 0",
                genCut, xMin, xMax, useVariableBinning);
  vGraphs.push_back(loose5);
  vLabels.push_back("p_{T} > 0");
  vColors.push_back(kGreen+2);

  TGraphAsymmErrors* medium5 = calculateEfficiency("genPhi", treePath, rootFileDirectory,
                 l1Cut + "&& pf_cPt > 2",
                 genCut, xMin, xMax, useVariableBinning);
  vGraphs.push_back(medium5);
  vLabels.push_back("p_{T} > 2");
  vColors.push_back(kBlue);

  TGraphAsymmErrors* tight5 = calculateEfficiency("genPhi", treePath, rootFileDirectory,
                l1Cut + "&& pf_cPt > 5",
                genCut, xMin, xMax, useVariableBinning);
  vGraphs.push_back(tight5);
  vLabels.push_back("p_{T} > 5");
  vColors.push_back(kBlack);

  TGraphAsymmErrors* all5 = calculateEfficiency("genPhi", treePath, rootFileDirectory,
              l1Cut + "&& pf_cPt > 10",
              genCut, xMin, xMax, useVariableBinning);
  vGraphs.push_back(all5);
  vLabels.push_back("p_{T} > 5");
  vColors.push_back(kRed);

  plotNEfficiencies(vGraphs, vLabels, vColors,
        "Gen #phi",
        "Phase 2 GCT",
        "efficiency_genPhi_barrel_PF_lowPt",
        outputDirectory);
  
}
/*********************************************************************/
