// Usage: root -l -b -q runValidationPlots.C
// Create plots of variables from a TTree, overlaying histograms from two cuts

#include "../baseCodeForPlots/comparisonPlots.C"

void runValidationPlots()
{
  // Load the macro
  gROOT->ProcessLine(".L ../baseCodeForPlots/comparisonPlots.C");
 
  TString treePath = "mutau_tree";
  TString inputDirectory  = "/eos/user/s/skkwan/hToAA/genSkims/Nov-03-2020-20h40m/SUSYVBFToHToAA_AToBB_AToTauTau_M-40GenSkim.root";
  TString outputDirectory = "plots/";

  TString sigCut = "";
  TString bkgCut = "";

  comparisonPlots("pt_1", sigCut, bkgCut, treePath, inputDirectory, outputDirectory, "pt_1", 60, 0, 120);
  comparisonPlots("eta_1", sigCut, bkgCut, treePath, inputDirectory, outputDirectory, "eta_1", 60, -3, 3);
  comparisonPlots("phi_1", sigCut, bkgCut, treePath, inputDirectory, outputDirectory, "phi_1", 60, -4, 4);

  comparisonPlots("pt_2", sigCut, bkgCut, treePath, inputDirectory, outputDirectory, "pt_2", 60, 0, 120);
  comparisonPlots("eta_2", sigCut, bkgCut, treePath, inputDirectory, outputDirectory, "eta_2", 60, -3, 3);
  comparisonPlots("phi_2", sigCut, bkgCut, treePath, inputDirectory, outputDirectory, "phi_2", 60, -4, 4);

  
}
