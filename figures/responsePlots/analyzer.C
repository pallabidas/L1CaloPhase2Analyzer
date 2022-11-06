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


   TFile hfile("hfile.root","RECREATE");
   TH1F *reso = new TH1F("reso", "reso", 40, 0., 2.);
   TH1F *reso5 = new TH1F("reso5", "reso5", 40, 0., 2.);
   TH1F *reso10 = new TH1F("reso10", "reso10", 40, 0., 2.);
   TH1F *reso15 = new TH1F("reso15", "reso15", 40, 0., 2.);
   TH1F *reso20 = new TH1F("reso20", "reso20", 40, 0., 2.);
   TH1F *reso25 = new TH1F("reso25", "reso25", 40, 0., 2.);
   TH1F *reso30 = new TH1F("reso30", "reso30", 40, 0., 2.);
   TH1F *reso35 = new TH1F("reso35", "reso35", 40, 0., 2.);
   TH1F *reso40 = new TH1F("reso40", "reso40", 40, 0., 2.);
   TH1F *reso45 = new TH1F("reso45", "reso45", 40, 0., 2.);
   TH1F *reso50 = new TH1F("reso50", "reso50", 40, 0., 2.);
   TH1F *reso55 = new TH1F("reso55", "reso55", 40, 0., 2.);
   TH1F *reso60 = new TH1F("reso60", "reso60", 40, 0., 2.);
   TH1F *reso70 = new TH1F("reso70", "reso70", 40, 0., 2.);
   TH1F *reso80 = new TH1F("reso80", "reso80", 40, 0., 2.);
   TH1F *reso90 = new TH1F("reso90", "reso90", 40, 0., 2.);
   TH1F *reso100 = new TH1F("reso100", "reso100", 40, 0., 2.);
   TH1F *reso120 = new TH1F("reso120", "reso120", 40, 0., 2.);
   TH1F *reso140 = new TH1F("reso140", "reso140", 40, 0., 2.);
   TH1F *reso160 = new TH1F("reso160", "reso160", 40, 0., 2.);

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;

   
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      //if(pf_cPt == 0) continue;
      reso->Fill(pf_cPt/genPt);
      if(genPt > 5. && genPt <= 10.) reso5->Fill(pf_cPt/genPt);
      if(genPt > 10. && genPt <= 15.) reso10->Fill(pf_cPt/genPt);
      if(genPt > 15. && genPt <= 20.) reso15->Fill(pf_cPt/genPt);
      if(genPt > 20. && genPt <= 25.) reso20->Fill(pf_cPt/genPt);
      if(genPt > 25. && genPt <= 30.) reso25->Fill(pf_cPt/genPt);
      if(genPt > 30. && genPt <= 35.) reso30->Fill(pf_cPt/genPt);
      if(genPt > 35. && genPt <= 40.) reso35->Fill(pf_cPt/genPt);
      if(genPt > 40. && genPt <= 45.) reso40->Fill(pf_cPt/genPt);
      if(genPt > 45. && genPt <= 50.) reso45->Fill(pf_cPt/genPt);
      if(genPt > 50. && genPt <= 55.) reso50->Fill(pf_cPt/genPt);
      if(genPt > 55. && genPt <= 60.) reso55->Fill(pf_cPt/genPt);
      if(genPt > 60. && genPt <= 70.) reso60->Fill(pf_cPt/genPt);
      if(genPt > 70. && genPt <= 80.) reso70->Fill(pf_cPt/genPt);
      if(genPt > 80. && genPt <= 90.) reso80->Fill(pf_cPt/genPt);
      if(genPt > 90. && genPt <= 100.) reso90->Fill(pf_cPt/genPt);
      if(genPt > 100. && genPt <= 120.) reso100->Fill(pf_cPt/genPt);
      if(genPt > 120. && genPt <= 140.) reso120->Fill(pf_cPt/genPt);
      if(genPt > 140. && genPt <= 160.) reso140->Fill(pf_cPt/genPt);
      if(genPt > 160.) reso160->Fill(pf_cPt/genPt);
   }

   hfile.Write();
   hfile.Close();
}
