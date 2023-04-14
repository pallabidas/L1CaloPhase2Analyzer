# setup.sh: get ROOT 5.34/24 (needs ROOT 5)

source /cvmfs/sft.cern.ch/lcg/views/LCG_71/x86_64-slc6-gcc48-opt/setup.sh 

cat > runDisplay.C << EOF
#include "TROOT.h"

void runDisplay() {

  gROOT->ProcessLine(".L plotEventDisplayPhaseIICaloJets.C++");

  gROOT->ProcessLine(".x plotEventDisplayPhaseIICaloJets.C(0)");


}

EOF

root -l -b -q runDisplay.C
