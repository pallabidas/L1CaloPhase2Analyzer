// Usage: root -l -b -q runDistributionPlots.C
// Create plots of variables from a TTree

#include "../baseCodeForPlots/singleDistribution.C"

void runDistributionPlots(TString sampleName, TString legend, TString inputDirectory)
{
  // Load the macro
  //  gROOT->ProcessLine(".L ../baseCodeForPlots/comparisonPlots.C");
 
  std::cout << sampleName << " " << legend << " " << inputDirectory << std::endl;
  TString treePath = "l1NtupleProducer/efficiencyTree";
  //  TString inputDirectory  = "/eos/user/s/skkwan/signalNanoAOD/2018/SUSYVBFToHToAA_AToBB_AToTauTau_M-40/SUSYVBFToHToAA_AToBB_AToTauTau_M-40_BATCH_1_NANO.root";
  TString outputDirectory = "/Users/stephaniekwan/Dropbox/Princeton_G4/Phase2RCT/analyzer/figures" + sampleName + "/"; 
  gSystem->Exec("mkdir -p " + outputDirectory);

  TString cut = "(genPt > 0)";  

  int nBins = 35;

  float ymaxDummyValue = -99;
  float ymaxPtDiff = 0.8;
  float ymaxDeltaR = 0.8;
  float ymaxPt = 0.14;
  float ymaxPhi = 0.05;

  // Gen quantities
  singleDistributionPlots("genPt",  "genPt",  cut, legend, treePath, inputDirectory, outputDirectory, "Gen Ele p_{T}", nBins, 0, 100,   0.06);
  singleDistributionPlots("genPhi", "genPhi", cut, legend, treePath, inputDirectory, outputDirectory, "Gen Ele #phi", nBins, -3.2, 3.2, ymaxDummyValue);
  singleDistributionPlots("genEta", "genEta", cut, legend, treePath, inputDirectory, outputDirectory, "Gen Ele #eta", nBins, -1.5, 1.5, ymaxDummyValue);
  

  // RCT
  singleDistributionPlots("rct_deltaR", "rct_deltaR", cut, legend, treePath, inputDirectory, outputDirectory, "RCT #DeltaR (L1, Gen)", nBins, 0, 0.5, ymaxDeltaR);
  singleDistributionPlots("rct_pT_fractional_diff", "(rct_cPt - genPt)/(genPt)", cut, legend, treePath, inputDirectory, outputDirectory,
                          "RCT (p_{T}^{reco} - p_{T}^{gen})/p_{T}^{gen}", nBins, -1.5, 1.5, ymaxPtDiff);

  singleDistributionPlots("rct_cPt",  "rct_cPt",  cut, legend, treePath, inputDirectory, outputDirectory, "RCT Cluster p_{T}", nBins, 0, 100, ymaxPt);
  singleDistributionPlots("rct_cPhi", "rct_cPhi", cut, legend, treePath, inputDirectory, outputDirectory, "RCT Cluster #phi", nBins, -3.2, 3.2, ymaxPhi);
  singleDistributionPlots("rct_cEta", "rct_cEta", cut, legend, treePath, inputDirectory, outputDirectory, "RCT Cluster #eta", nBins, -1.5, 1.5, ymaxDummyValue);


  // GCT
  singleDistributionPlots("gct_deltaR", "gct_deltaR", cut, legend, treePath, inputDirectory, outputDirectory, "GCT #DeltaR (L1, Gen)", nBins, 0, 0.5, ymaxDeltaR);
  singleDistributionPlots("gct_pT_fractional_diff", "(gct_cPt - genPt)/(genPt)", cut, legend, treePath, inputDirectory, outputDirectory,
                          "GCT (p_{T}^{reco} - p_{T}^{gen})/p_{T}^{gen}", nBins, -1.5, 1.5, ymaxPtDiff);

  singleDistributionPlots("gct_cPt",  "gct_cPt",  cut, legend, treePath, inputDirectory, outputDirectory, "GCT Cluster p_{T}", nBins, 0, 100, ymaxPt);
  singleDistributionPlots("gct_cPhi", "gct_cPhi", cut, legend, treePath, inputDirectory, outputDirectory, "GCT Cluster #phi", nBins, -3.2, 3.2, ymaxPhi);
  singleDistributionPlots("gct_cEta", "gct_cEta", cut, legend, treePath, inputDirectory, outputDirectory, "GCT Cluster #eta", nBins, -1.5, 1.5, ymaxDummyValue);


}
