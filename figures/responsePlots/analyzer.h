//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Nov  5 19:51:39 2022 by ROOT version 6.22/09
// from TTree efficiencyTree/Efficiency Tree
// found on file: analyzer_singlepion.root
//////////////////////////////////////////////////////////

#ifndef analyzer_h
#define analyzer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class analyzer {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   vector<TLorentzVector> *rctClusters;
   vector<TLorentzVector> *rctTowers;
   vector<TLorentzVector> *hcalTPGs;
   vector<TLorentzVector> *ecalTPGs;
   vector<TLorentzVector> *gctTowers;
   vector<TLorentzVector> *caloPFClusters;
   Int_t           run;
   Int_t           lumi;
   Int_t           event;
   Int_t           nvtx;
   Double_t        genPt;
   Double_t        genEta;
   Double_t        genPhi;
   Double_t        rct_cPt;
   Double_t        rct_cEta;
   Double_t        rct_cPhi;
   Double_t        rct_deltaR;
   Double_t        rct_et2x5;
   Double_t        rct_et5x5;
   Double_t        gct_cPt;
   Double_t        gct_cEta;
   Double_t        gct_cPhi;
   Double_t        gct_deltaR;
   Double_t        gct_et2x5;
   Double_t        gct_et5x5;
   Double_t        gct_iso;
   Int_t           gct_is_ss;
   Int_t           gct_is_looseTkss;
   Int_t           gct_is_iso;
   Int_t           gct_is_looseTkiso;
   Double_t        pf_cPt;
   Double_t        pf_cEta;
   Double_t        pf_cPhi;
   Double_t        pf_deltaR;

   // List of branches
   TBranch        *b_rctClusters;   //!
   TBranch        *b_rctTowers;   //!
   TBranch        *b_hcalTPGs;   //!
   TBranch        *b_ecalTPGs;   //!
   TBranch        *b_gctTowers;   //!
   TBranch        *b_caloPFClusters;   //!
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_event;   //!
   TBranch        *b_nvtx;   //!
   TBranch        *b_genPt;   //!
   TBranch        *b_genEta;   //!
   TBranch        *b_genPhi;   //!
   TBranch        *b_rct_cPt;   //!
   TBranch        *b_rct_cEta;   //!
   TBranch        *b_rct_cPhi;   //!
   TBranch        *b_rct_deltaR;   //!
   TBranch        *b_rct_et2x5;   //!
   TBranch        *b_rct_et5x5;   //!
   TBranch        *b_gct_cPt;   //!
   TBranch        *b_gct_cEta;   //!
   TBranch        *b_gct_cPhi;   //!
   TBranch        *b_gct_deltaR;   //!
   TBranch        *b_gct_et2x5;   //!
   TBranch        *b_gct_et5x5;   //!
   TBranch        *b_gct_iso;   //!
   TBranch        *b_gct_is_ss;   //!
   TBranch        *b_gct_is_looseTkss;   //!
   TBranch        *b_gct_is_iso;   //!
   TBranch        *b_gct_is_looseTkiso;   //!
   TBranch        *b_pf_cPt;   //!
   TBranch        *b_pf_cEta;   //!
   TBranch        *b_pf_cPhi;   //!
   TBranch        *b_pf_deltaR;   //!

   analyzer(TTree *tree=0);
   virtual ~analyzer();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef analyzer_cxx
analyzer::analyzer(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("analyzer_singlepion.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("analyzer_singlepion.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("analyzer_singlepion.root:/l1NtupleProducer");
      dir->GetObject("efficiencyTree",tree);

   }
   Init(tree);
}

analyzer::~analyzer()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t analyzer::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t analyzer::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void analyzer::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   rctClusters = 0;
   rctTowers = 0;
   hcalTPGs = 0;
   ecalTPGs = 0;
   gctTowers = 0;
   caloPFClusters = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("rctClusters", &rctClusters, &b_rctClusters);
   fChain->SetBranchAddress("rctTowers", &rctTowers, &b_rctTowers);
   fChain->SetBranchAddress("hcalTPGs", &hcalTPGs, &b_hcalTPGs);
   fChain->SetBranchAddress("ecalTPGs", &ecalTPGs, &b_ecalTPGs);
   fChain->SetBranchAddress("gctTowers", &gctTowers, &b_gctTowers);
   fChain->SetBranchAddress("caloPFClusters", &caloPFClusters, &b_caloPFClusters);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("nvtx", &nvtx, &b_nvtx);
   fChain->SetBranchAddress("genPt", &genPt, &b_genPt);
   fChain->SetBranchAddress("genEta", &genEta, &b_genEta);
   fChain->SetBranchAddress("genPhi", &genPhi, &b_genPhi);
   fChain->SetBranchAddress("rct_cPt", &rct_cPt, &b_rct_cPt);
   fChain->SetBranchAddress("rct_cEta", &rct_cEta, &b_rct_cEta);
   fChain->SetBranchAddress("rct_cPhi", &rct_cPhi, &b_rct_cPhi);
   fChain->SetBranchAddress("rct_deltaR", &rct_deltaR, &b_rct_deltaR);
   fChain->SetBranchAddress("rct_et2x5", &rct_et2x5, &b_rct_et2x5);
   fChain->SetBranchAddress("rct_et5x5", &rct_et5x5, &b_rct_et5x5);
   fChain->SetBranchAddress("gct_cPt", &gct_cPt, &b_gct_cPt);
   fChain->SetBranchAddress("gct_cEta", &gct_cEta, &b_gct_cEta);
   fChain->SetBranchAddress("gct_cPhi", &gct_cPhi, &b_gct_cPhi);
   fChain->SetBranchAddress("gct_deltaR", &gct_deltaR, &b_gct_deltaR);
   fChain->SetBranchAddress("gct_et2x5", &gct_et2x5, &b_gct_et2x5);
   fChain->SetBranchAddress("gct_et5x5", &gct_et5x5, &b_gct_et5x5);
   fChain->SetBranchAddress("gct_iso", &gct_iso, &b_gct_iso);
   fChain->SetBranchAddress("gct_is_ss", &gct_is_ss, &b_gct_is_ss);
   fChain->SetBranchAddress("gct_is_looseTkss", &gct_is_looseTkss, &b_gct_is_looseTkss);
   fChain->SetBranchAddress("gct_is_iso", &gct_is_iso, &b_gct_is_iso);
   fChain->SetBranchAddress("gct_is_looseTkiso", &gct_is_looseTkiso, &b_gct_is_looseTkiso);
   fChain->SetBranchAddress("pf_cPt", &pf_cPt, &b_pf_cPt);
   fChain->SetBranchAddress("pf_cEta", &pf_cEta, &b_pf_cEta);
   fChain->SetBranchAddress("pf_cPhi", &pf_cPhi, &b_pf_cPhi);
   fChain->SetBranchAddress("pf_deltaR", &pf_deltaR, &b_pf_deltaR);
   Notify();
}

Bool_t analyzer::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void analyzer::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t analyzer::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef analyzer_cxx
