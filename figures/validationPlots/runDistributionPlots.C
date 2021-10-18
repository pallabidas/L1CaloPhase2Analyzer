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

  int nBins = 50;
  singleDistributionPlots("deltaR", cut, legend, treePath, inputDirectory, outputDirectory, "#DeltaR (L1, Gen)", nBins, 0, 0.5);


  // singleDistributionPlots("pt_2", cut, legend, treePath, inputDirectory, outputDirectory, "#tau_{h} p_{T}", nBins, 15, 150);
  // singleDistributionPlots("eta_2", cut, legend, treePath, inputDirectory, outputDirectory, "#tau_{h} #eta", nBins, -3, 3);
  // singleDistributionPlots("phi_2", cut, legend, treePath, inputDirectory, outputDirectory, "#tau_{h} #phi", nBins, -4, 4);

  

  
}
