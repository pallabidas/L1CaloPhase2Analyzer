// Usage: root -l -b -q runDistributionPlots.C
// Create plots of variables from a TTree

#include "../baseCodeForPlots/singleDistribution.C"
#include "../baseCodeForPlots/scatterPlot.cpp"

void makeDistributionPlots(TString sampleName, TString legend, TString inputDirectory)
{
  // Load the macro
  //  gROOT->ProcessLine(".L ../baseCodeForPlots/comparisonPlots.C");
 
  std::cout << sampleName << " " << legend << " " << inputDirectory << std::endl;
  TString treePath = "l1NtupleProducer/efficiencyTree";
  //  TString inputDirectory  = "/eos/user/s/skkwan/signalNanoAOD/2018/SUSYVBFToHToAA_AToBB_AToTauTau_M-40/SUSYVBFToHToAA_AToBB_AToTauTau_M-40_BATCH_1_NANO.root";
  TString outputDirectory = "/Users/stephaniekwan/Dropbox/Princeton_G4/Phase2RCT/analyzer/figures" + sampleName + "/"; 
  gSystem->Exec("mkdir -p " + outputDirectory);

  TString cut = "(genPt > 0) && (abs(genEta) < 1.305)";  
  TString defaultBonusDescriptor = "";
  TString tdrCutDescriptor = "L1 p_{T} > 15";

  int nBins = 35;

  float ymaxDummyValue = -99;
  float ymaxPtDiff = 0.8;
  float ymaxDeltaR = 0.3;
  float ymaxPt = 0.06;
  float ymaxPhi = 0.05;

  // // Troubleshooting:
  // // 
  // TString poorMatchCut = "(gct_cPt - genPt)/(genPt) < -0.6";  
  // TString bonusDescriptorPoorMatchCut = "#scale[0.6]{p_{T} diff < -0.6}";
  // scatterPlot("pTdiff_lt_neg0p6_scatter_gct_cPt_vs_genPt", "gct_cPt:genPt", poorMatchCut, legend, treePath, inputDirectory, outputDirectory, 
  //             "Gen Ele p_{T}", "GCT Cluster p_{T}",  bonusDescriptorPoorMatchCut, nBins, 0, 100, nBins, 0, 100);
  // singleDistributionPlots("pTdiff_lt_neg0p6_genPt",  "genPt",  poorMatchCut, legend, treePath, inputDirectory, outputDirectory, "Gen Ele p_{T}", bonusDescriptorPoorMatchCut, nBins, 0, 100, 0.10);
  // singleDistributionPlots("pTdiff_lt_neg0p6_genPhi", "genPhi", poorMatchCut, legend, treePath, inputDirectory, outputDirectory, "Gen Ele #phi", bonusDescriptorPoorMatchCut, nBins, -3.2, 3.2, ymaxDummyValue);
  // singleDistributionPlots("pTdiff_lt_neg0p6_genEta", "genEta", poorMatchCut, legend, treePath, inputDirectory, outputDirectory, "Gen Ele #eta", bonusDescriptorPoorMatchCut, nBins, -1.5, 1.5, ymaxDummyValue);
  // singleDistributionPlots("pTdiff_lt_neg0p6_gct_deltaR", "gct_deltaR", poorMatchCut, legend, treePath, inputDirectory, outputDirectory, "GCT #DeltaR (L1, Gen)", bonusDescriptorPoorMatchCut, nBins, 0, 0.5, ymaxDummyValue);
  // singleDistributionPlots("pTdiff_lt_neg0p6_gct_pT_fractional_diff", "(gct_cPt - genPt)/(genPt)", poorMatchCut, legend, treePath, inputDirectory, outputDirectory,
  //                          "GCT (p_{T}^{reco} - p_{T}^{gen})/p_{T}^{gen}", bonusDescriptorPoorMatchCut, nBins, -1.5, 1.5, ymaxPtDiff);


  // TString goodMatchCut = "(gct_cPt - genPt)/(genPt) > -0.6";  
  // TString bonusDescriptorGoodMatchCut = "#scale[0.6]{p_{T} diff > -0.6}";
  // scatterPlot("pTdiff_gr_neg0p6_scatter_gct_cPt_vs_genPt", "gct_cPt:genPt", goodMatchCut, legend, treePath, inputDirectory, outputDirectory, 
  //             "Gen Ele p_{T}", "GCT Cluster p_{T}",  bonusDescriptorGoodMatchCut, nBins, 0, 100, nBins, 0, 100);
  // singleDistributionPlots("pTdiff_gr_neg0p6_genPt",  "genPt",  goodMatchCut, legend, treePath, inputDirectory, outputDirectory, "Gen Ele p_{T}", bonusDescriptorGoodMatchCut, nBins, 0, 100, 0.10);
  // singleDistributionPlots("pTdiff_gr_neg0p6_genPhi", "genPhi", goodMatchCut, legend, treePath, inputDirectory, outputDirectory, "Gen Ele #phi", bonusDescriptorGoodMatchCut, nBins, -3.2, 3.2, ymaxDummyValue);
  // singleDistributionPlots("pTdiff_gr_neg0p6_genEta", "genEta", goodMatchCut, legend, treePath, inputDirectory, outputDirectory, "Gen Ele #eta", bonusDescriptorGoodMatchCut, nBins, -1.5, 1.5, ymaxDummyValue);
  // singleDistributionPlots("pTdiff_gr_neg0p6_gct_deltaR", "gct_deltaR", goodMatchCut, legend, treePath, inputDirectory, outputDirectory, "GCT #DeltaR (L1, Gen)", bonusDescriptorGoodMatchCut, nBins, 0, 0.5, ymaxDummyValue);
  // singleDistributionPlots("pTdiff_gr_neg0p6_gct_pT_fractional_diff", "(gct_cPt - genPt)/(genPt)", goodMatchCut, legend, treePath, inputDirectory, outputDirectory,
  //                          "GCT (p_{T}^{reco} - p_{T}^{gen})/p_{T}^{gen}", bonusDescriptorGoodMatchCut, nBins, -1.5, 1.5, ymaxPtDiff);


  // Gen quantities
  singleDistributionPlots("genPt",  "genPt",  cut, legend, treePath, inputDirectory, outputDirectory, "Gen Ele p_{T}", defaultBonusDescriptor, nBins, 0, 100,   0.06);
  singleDistributionPlots("genPhi", "genPhi", cut, legend, treePath, inputDirectory, outputDirectory, "Gen Ele #phi", defaultBonusDescriptor, nBins, -3.2, 3.2, ymaxDummyValue);
  singleDistributionPlots("genEta", "genEta", cut, legend, treePath, inputDirectory, outputDirectory, "Gen Ele #eta", defaultBonusDescriptor, nBins, -1.5, 1.5, ymaxDummyValue);
  

  // RCT
  // Note that the next two plots have an additional pT cut (reflected in the caption)
  singleDistributionPlots("rct_deltaR", "rct_deltaR", cut + "&& (rct_cPt > 15)", legend, treePath, inputDirectory, outputDirectory, "RCT #DeltaR (L1, Gen)", tdrCutDescriptor, nBins, 0, 0.10, ymaxDeltaR);
  singleDistributionPlots("rct_pT_fractional_diff", "(rct_cPt - genPt)/(genPt)", cut + "&& (rct_cPt > 15)", legend, treePath, inputDirectory, outputDirectory,
                          "RCT (p_{T}^{reco} - p_{T}^{gen})/p_{T}^{gen}", tdrCutDescriptor, nBins, -1.2, 1.2, ymaxPtDiff);

  singleDistributionPlots("rct_cPt",  "rct_cPt",  cut, legend, treePath, inputDirectory, outputDirectory, "RCT Cluster p_{T}", defaultBonusDescriptor, nBins, 0, 100, ymaxPt);
  singleDistributionPlots("rct_cPhi", "rct_cPhi", cut, legend, treePath, inputDirectory, outputDirectory, "RCT Cluster #phi", defaultBonusDescriptor, nBins, -3.2, 3.2, ymaxPhi);
  singleDistributionPlots("rct_cEta", "rct_cEta", cut, legend, treePath, inputDirectory, outputDirectory, "RCT Cluster #eta", defaultBonusDescriptor, nBins, -1.5, 1.5, ymaxDummyValue);


  // GCT
  // Note that the next two plots have an additional pT cut (reflected in the caption)
  singleDistributionPlots("gct_deltaR", "gct_deltaR", cut + "&& (gct_cPt > 15)", legend, treePath, inputDirectory, outputDirectory, "GCT #DeltaR (L1, Gen)", tdrCutDescriptor, nBins, 0, 0.10, ymaxDeltaR);
  singleDistributionPlots("gct_pT_fractional_diff", "(gct_cPt - genPt)/(genPt)", cut + "&& (gct_cPt > 15)", legend, treePath, inputDirectory, outputDirectory,
                          "GCT (p_{T}^{reco} - p_{T}^{gen})/p_{T}^{gen}", tdrCutDescriptor, nBins, -1.2, 1.2, ymaxPtDiff);

  singleDistributionPlots("gct_cPt",  "gct_cPt",  cut, legend, treePath, inputDirectory, outputDirectory, "GCT Cluster p_{T}", defaultBonusDescriptor, nBins, 0, 100, ymaxPt);
  singleDistributionPlots("gct_cPhi", "gct_cPhi", cut, legend, treePath, inputDirectory, outputDirectory, "GCT Cluster #phi", defaultBonusDescriptor, nBins, -3.2, 3.2, ymaxPhi);
  singleDistributionPlots("gct_cEta", "gct_cEta", cut, legend, treePath, inputDirectory, outputDirectory, "GCT Cluster #eta", defaultBonusDescriptor, nBins, -1.5, 1.5, ymaxDummyValue);


}
