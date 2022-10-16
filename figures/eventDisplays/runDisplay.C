#include "TROOT.h"

void runDisplay() {

  gROOT->ProcessLine(".L plotEventDisplayPhaseIIPFclusters.C++");

  gROOT->ProcessLine(".x plotEventDisplayPhaseIIPFclusters.C(0)");


}

