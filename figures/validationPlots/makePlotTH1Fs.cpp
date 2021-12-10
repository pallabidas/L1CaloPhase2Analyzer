// Usage: root -l -b -q makePlotTH1F.cpp

// Customize this file

#include "../baseCodeForPlots/plotTH1F.cpp"

void makePlotTH1Fs()
{
  

/* Read in and plot an EXISTING TH1F:
    th1fName :  the name of the TH1F inside the file (usually the name of the variable).
    xLabel   :  the x-axis label on the plot.
    paveText :  text to go in a TPave box. 
    isAU     : if True, normalize the plot so the area under the curve is 1, and set y-axis title to "A.U.".
               if False, do not normalize and set y-axis title to "Counts".
    histPath :  the path to the TH1F inside the file (usually the name of the folder in the ROOT file).
    inputDirectory : the path to the ROOT file (~/myRepo/input.root).
    outputDirectory: the directory where the output plots will be saved.
    */
 
  TString treePath = "l1NtupleProducer/";
  TString inputDirectory  = "/Users/stephaniekwan/Dropbox/Princeton_G4/Phase2RCT/analyzer/analyzer-rates.root";
  TString outputDirectory = "/Users/stephaniekwan/Dropbox/Princeton_G4/Phase2RCT/analyzer/ratesPlots/";

  TString label = "MinBias 200 PU, no cuts";
  bool isAU = false;
  bool useLogy = false;
  plotTH1F("l1eg_pt", "L1 cluster p_{T}", label, isAU, treePath, inputDirectory, outputDirectory, true);
}
