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

  TString cut = "";

  int nBins = 35;
  singleDistributionPlots("deltaR", "deltaR", cut, legend, treePath, inputDirectory, outputDirectory, "#DeltaR (L1, Gen)", nBins, 0, 0.5);
  singleDistributionPlots("pT_fractional_diff", "(cPt - genPt)/(genPt)", cut, legend, treePath, inputDirectory, outputDirectory,
                          "(p_{T}^{reco} - p_{T}^{gen})/p_{T}^{gen}", nBins, -1.5, 1.5);

  singleDistributionPlots("cPt",  "cPt",  cut, legend, treePath, inputDirectory, outputDirectory, "Cluster p_{T}", nBins, 0, 100);
  singleDistributionPlots("cPhi", "cPhi", cut, legend, treePath, inputDirectory, outputDirectory, "Cluster #phi", nBins, -3.2, 3.2);
  singleDistributionPlots("cEta", "cEta", cut, legend, treePath, inputDirectory, outputDirectory, "Cluster #eta", nBins, -1.5, 1.5);

  singleDistributionPlots("genPt",  "genPt",  cut, legend, treePath, inputDirectory, outputDirectory, "Gen EG p_{T}", nBins, 0, 100);
  singleDistributionPlots("genPhi", "genPhi", cut, legend, treePath, inputDirectory, outputDirectory, "Gen EG #phi", nBins, -3.2, 3.2);
  singleDistributionPlots("genEta", "genEta", cut, legend, treePath, inputDirectory, outputDirectory, "Gen EG #eta", nBins, -1.5, 1.5);
  

  
}
