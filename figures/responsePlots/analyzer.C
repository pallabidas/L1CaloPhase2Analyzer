#define analyzer_cxx
#include "analyzer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void analyzer::Loop()
{
//   In a ROOT session, you can do:
//      root> .L analyzer.C
//      root> analyzer t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;


   TFile hfile("hfile_forCalib.root","RECREATE");
   TH1F *reso = new TH1F("reso", "reso", 40, 0., 200.);
   TH1F *reso5 = new TH1F("reso5", "reso5", 20, 0., 30.);
   TH1F *reso10 = new TH1F("reso10", "reso10", 20, 0., 40.);
   TH1F *reso15 = new TH1F("reso15", "reso15", 20, 0., 60.);
   TH1F *reso20 = new TH1F("reso20", "reso20", 20, 0., 60.);
   TH1F *reso25 = new TH1F("reso25", "reso25", 20, 5., 70.);
   TH1F *reso30 = new TH1F("reso30", "reso30", 20, 10., 80.);
   TH1F *reso35 = new TH1F("reso35", "reso35", 20, 15., 100.);
   TH1F *reso40 = new TH1F("reso40", "reso40", 20, 20., 110.);
   TH1F *reso45 = new TH1F("reso45", "reso45", 20, 25., 120.);
   TH1F *reso50 = new TH1F("reso50", "reso50", 20, 30., 120.);
   TH1F *reso55 = new TH1F("reso55", "reso55", 20, 35., 130.);
   TH1F *reso60 = new TH1F("reso60", "reso60", 20, 40., 150.);
   TH1F *reso70 = new TH1F("reso70", "reso70", 20, 50., 160.);
   TH1F *reso80 = new TH1F("reso80", "reso80", 20, 60., 170.);
   TH1F *reso90 = new TH1F("reso90", "reso90", 20, 70., 200.);
   TH1F *reso100 = new TH1F("reso100", "reso100", 20, 80., 200.);
   TH1F *reso120 = new TH1F("reso120", "reso120", 20, 100., 200.);
   TH1F *reso140 = new TH1F("reso140", "reso140", 20, 120., 200.);
   TH1F *reso160 = new TH1F("reso160", "reso160", 20, 140., 200.);

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      //if(pf_cPt == 0) continue;
      reso->Fill(pf_cPt);
      if(pf_cPt > 5. && pf_cPt <= 10.) reso5->Fill(genPt);
      if(pf_cPt > 10. && pf_cPt <= 15.) reso10->Fill(genPt);
      if(pf_cPt > 15. && pf_cPt <= 20.) reso15->Fill(genPt);
      if(pf_cPt > 20. && pf_cPt <= 25.) reso20->Fill(genPt);
      if(pf_cPt > 25. && pf_cPt <= 30.) reso25->Fill(genPt);
      if(pf_cPt > 30. && pf_cPt <= 35.) reso30->Fill(genPt);
      if(pf_cPt > 35. && pf_cPt <= 40.) reso35->Fill(genPt);
      if(pf_cPt > 40. && pf_cPt <= 45.) reso40->Fill(genPt);
      if(pf_cPt > 45. && pf_cPt <= 50.) reso45->Fill(genPt);
      if(pf_cPt > 50. && pf_cPt <= 55.) reso50->Fill(genPt);
      if(pf_cPt > 55. && pf_cPt <= 60.) reso55->Fill(genPt);
      if(pf_cPt > 60. && pf_cPt <= 70.) reso60->Fill(genPt);
      if(pf_cPt > 70. && pf_cPt <= 80.) reso70->Fill(genPt);
      if(pf_cPt > 80. && pf_cPt <= 90.) reso80->Fill(genPt);
      if(pf_cPt > 90. && pf_cPt <= 100.) reso90->Fill(genPt);
      if(pf_cPt > 100. && pf_cPt <= 120.) reso100->Fill(genPt);
      if(pf_cPt > 120. && pf_cPt <= 140.) reso120->Fill(genPt);
      if(pf_cPt > 140. && pf_cPt <= 160.) reso140->Fill(genPt);
      if(pf_cPt > 160.) reso160->Fill(genPt);
   }

   hfile.Write();
   hfile.Close();
}
