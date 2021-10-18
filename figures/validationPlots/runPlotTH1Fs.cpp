// Usage: root -l -b -q runPlotTH1F.cpp

// Customize this file

#include "../baseCodeForPlots/plotTH1F.cpp"

void runPlotTH1Fs()
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
 
  TString treePath = "mutau_slimmed/";
  TString inputDirectory  = "~/Dropbox/Princeton G3/Pre-Thesis/hists.root";
  TString outputDirectory = "myPlots/";

  TString label = "(MC) ggH #rightarrow aa #rightarrow bb#tau#tau, #mu#tau_{H} final state";
  TString llabel = "#mu#tau_{H}";
  bool isAU = false;


  plotTH1F("muPt", "Muon p_{T}", label, isAU, treePath, inputDirectory, outputDirectory);
  plotTH1F("muEta", "Muon #eta", label, isAU, treePath, inputDirectory, outputDirectory);
  plotTH1F("tauPt", "#tau_{H} p_{t}", label, isAU, treePath, inputDirectory, outputDirectory);
  plotTH1F("tauEta", "#tau_{H} #eta", label, isAU, treePath, inputDirectory, outputDirectory);
  
}
