// Usage: root -l -b -q runDistributionPlots.C
// Create plots of variables from a TTree

#include "../baseCodeForPlots/singleDistribution.C"
#include "../baseCodeForPlots/scatterPlot.cpp"

// sampleName: name of folder                      (e.g. "RelValElectronGunPt2To100_NoPU")
// legend: description at top of plot, LaTeX OK    (e.g. "RelVal Electron Gun Pt 2 to 100, No PU"")
// inputDirectory: path to input root file         (e.g. "analyzer.root")
// cut: cut applied to all plots                   (e.g. "pt > 0")
// descriptor: first descriptin in plot, LaTeX OK  (e.g. "L1 p_{T} > 0" )
// descriptor, bonusDescriptor: second line of description in plot, LaTeX OK 
// SF: scale factor for cluster pT (specific to this project)
void makeDistributionPlots(TString sampleName, TString legend, TString inputDirectory,
                           TString cut = "", 
                           TString descriptor = "",
                           TString bonusDescriptor = "",
                           TString SF = "1.0")
{
  // Load the macro
  //  gROOT->ProcessLine(".L ../baseCodeForPlots/comparisonPlots.C");
 
  std::cout << sampleName << " " << legend << " " << inputDirectory << std::endl;
  TString treePath = "l1NtupleProducer/efficiencyTree";
  //  TString inputDirectory  = "/eos/user/s/skkwan/signalNanoAOD/2018/SUSYVBFToHToAA_AToBB_AToTauTau_M-40/SUSYVBFToHToAA_AToBB_AToTauTau_M-40_BATCH_1_NANO.root";
  TString outputDirectory = "/Users/stephaniekwan/Dropbox/Princeton_G4/Phase2RCT/analyzer/studies/" + sampleName + "/"; 
  gSystem->Exec("mkdir -p " + outputDirectory);

  // TString cut = "(gct_cPt > 20) && (gct_cPt < 75) && (genPt > 30) && (genPt < 70)";    // changed from genPt > 0
  //  TString descriptor, bonusDescriptor = "L1 p_{T}/0.985, 20 < L1 p_{T} < 75, 30 < Gen p_{T} < 70";

  int nBins = 51;

  float ymaxDummyValue = -99;
  float ymaxPtDiff = 0.35;
  float ymaxDeltaR = 0.1;
  float ymaxPt = 0.06;
  float ymaxPhi = 0.05;

  descriptor = ("|#eta^{Gen}| < 1.4841, " + descriptor);  

  // GCT
  // Note that the next two plots have an additional pT cut (reflected in the caption)
  // singleDistributionPlots("gct_deltaR", "gct_deltaR", cut + "&& (gct_cPt > 15)" + poorMatchCut, legend, treePath, inputDirectory, outputDirectory, "GCT #DeltaR (L1, Gen)", descriptor, bonusDescriptorPoorMatchCut, nBins, 0, 0.10, ymaxDeltaR);
  // singleDistributionPlots("gct_pT_fractional_diff", "(gct_cPt - genPt)/(genPt)", cut + "&& (gct_cPt > 15)" + poorMatchCut, legend, treePath, inputDirectory, outputDirectory,
  //                         "GCT (p_{T}^{reco} - p_{T}^{gen})/p_{T}^{gen}", descriptor, bonusDescriptorPoorMatchCut, nBins, -0.5, 0.5, ymaxPtDiff);

  // singleDistributionPlots("gct_cPt",  "gct_cPt",  cut + "&& (gct_cPt > 15)" + poorMatchCut, legend, treePath, inputDirectory, outputDirectory, "GCT Cluster p_{T}", descriptor, bonusDescriptorPoorMatchCut, nBins, 0, 100, 0.2);
  // singleDistributionPlots("gct_cPhi", "gct_cPhi", cut + "&& (gct_cPt > 15)" + poorMatchCut, legend, treePath, inputDirectory, outputDirectory, "GCT Cluster #phi", descriptor, bonusDescriptorPoorMatchCut, nBins, -3.2, 3.2, 0.1);
  // singleDistributionPlots("gct_cEta", "gct_cEta", cut + "&& (gct_cPt > 15)" + poorMatchCut, legend, treePath, inputDirectory, outputDirectory, "GCT Cluster #eta", descriptor, bonusDescriptorPoorMatchCut, nBins, -1.5, 1.5, ymaxDummyValue);

  // singleDistributionPlots("genPt",  "genPt",  cut + "&& (gct_cPt > 15)" + poorMatchCut, legend, treePath, inputDirectory, outputDirectory, "Gen Ele p_{T}", descriptor, bonusDescriptorPoorMatchCut, nBins, 0, 100,   0.2);
  // singleDistributionPlots("genPhi", "genPhi", cut + "&& (gct_cPt > 15)" + poorMatchCut, legend, treePath, inputDirectory, outputDirectory, "Gen Ele #phi", descriptor, bonusDescriptorPoorMatchCut, nBins, -3.2, 3.2, 0.1);
  // singleDistributionPlots("genEta", "genEta", cut + "&& (gct_cPt > 15)" + poorMatchCut, legend, treePath, inputDirectory, outputDirectory, "Gen Ele #eta", descriptor, bonusDescriptorPoorMatchCut, nBins, -1.5, 1.5, ymaxDummyValue);
  

  // Gen quantities
  singleDistributionPlots("genPt",  "genPt",  cut, legend, treePath, inputDirectory, outputDirectory, "Gen Ele p_{T}", descriptor, bonusDescriptor, nBins, 0, 100,   0.06);
  singleDistributionPlots("genPhi", "genPhi", cut, legend, treePath, inputDirectory, outputDirectory, "Gen Ele #phi", descriptor, bonusDescriptor, nBins, -3.2, 3.2, ymaxDummyValue);
  singleDistributionPlots("genEta", "genEta", cut, legend, treePath, inputDirectory, outputDirectory, "Gen Ele #eta", descriptor, bonusDescriptor, nBins, -1.5, 1.5, ymaxDummyValue);
  

  // RCT
  // Note that the next two plots have an additional pT cut (reflected in the caption)
  // singleDistributionPlots("rct_deltaR", "rct_deltaR", cut + "&& (rct_cPt > 15)", legend, treePath, inputDirectory, outputDirectory, "RCT #DeltaR (L1, Gen)", descriptor, bonusDescriptor, nBins, 0, 0.10, ymaxDeltaR);
  // singleDistributionPlots("rct_pT_fractional_diff", "(rct_cPt - genPt)/(genPt)", cut + "&& (rct_cPt > 15)", legend, treePath, inputDirectory, outputDirectory,
  //                        "RCT (p_{T}^{reco} - p_{T}^{gen})/p_{T}^{gen}", descriptor, bonusDescriptor, nBins, -0.5, 0.5, ymaxPtDiff);

  // singleDistributionPlots("rct_cPt",  "rct_cPt",  cut, legend, treePath, inputDirectory, outputDirectory, "RCT Cluster p_{T}", defaultdescriptor, bonusDescriptor, nBins, 0, 100, ymaxPt);
  // singleDistributionPlots("rct_cPhi", "rct_cPhi", cut, legend, treePath, inputDirectory, outputDirectory, "RCT Cluster #phi", defaultdescriptor, bonusDescriptor, nBins, -3.2, 3.2, ymaxPhi);
  // singleDistributionPlots("rct_cEta", "rct_cEta", cut, legend, treePath, inputDirectory, outputDirectory, "RCT Cluster #eta", defaultdescriptor, bonusDescriptor, nBins, -1.5, 1.5, ymaxDummyValue);


  // GCT
  singleDistributionPlots("gct_deltaR", "gct_deltaR", cut + "", legend, treePath, inputDirectory, outputDirectory, "GCT #DeltaR (L1, Gen)", descriptor, bonusDescriptor, nBins, 0, 0.05, ymaxDeltaR);
  singleDistributionPlots("gct_pT_resolution", "((gct_cPt/" + SF + ") - genPt)/(genPt)", cut + "", legend, treePath, inputDirectory, outputDirectory,
                         "GCT (p_{T}^{reco} - p_{T}^{gen})/p_{T}^{gen}", descriptor, bonusDescriptor, nBins, -0.5, 0.5, ymaxPtDiff);

  singleDistributionPlots("gct_cPt",  "(gct_cPt/" + SF + ")",  cut, legend, treePath, inputDirectory, outputDirectory, "GCT Cluster p_{T}", descriptor, bonusDescriptor, nBins, 0, 100, ymaxPt);
  singleDistributionPlots("gct_cPhi", "gct_cPhi", cut, legend, treePath, inputDirectory, outputDirectory, "GCT Cluster #phi", descriptor, bonusDescriptor, nBins, -3.2, 3.2, ymaxPhi);
  singleDistributionPlots("gct_cEta", "gct_cEta", cut, legend, treePath, inputDirectory, outputDirectory, "GCT Cluster #eta", descriptor, bonusDescriptor, nBins, -1.5, 1.5, ymaxDummyValue);

  singleDistributionPlots("gct_eta_diff", "genEta - gct_cEta", cut, legend, treePath, inputDirectory, outputDirectory, "Gen #eta - GCT Cluster #eta", descriptor, bonusDescriptor, nBins, -0.06, 0.06, ymaxDummyValue);
  singleDistributionPlots("gct_phi_diff", "genPhi - gct_cPhi", cut, legend, treePath, inputDirectory, outputDirectory, "Gen #phi - GCT Cluster #phi", descriptor, bonusDescriptor, nBins, -0.06, 0.06, ymaxDummyValue);
}
