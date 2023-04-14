#include "TROOT.h"

void runDisplay() {

  gROOT->ProcessLine(".L plotEventDisplayPhaseIICaloJets.C++");

  gROOT->ProcessLine(".x plotEventDisplayPhaseIICaloJets.C(0)");


}

