{
  TCanvas* Tcan = new TCanvas("Tcan","", 100, 20, 1000, 800);
  ifstream f1;
  f1.open("scale_pfcluster.txt");
  string s;
  float yvalue, yerr;
  Float_t bins[] = {5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100, 120, 150, 200};
  Int_t  binnum = 18;
  TH1F* h1 = new TH1F("h1","h1", binnum, bins);

  for(int i = 0; i < binnum; i++){
    f1>>s>>yvalue>>yerr;
    h1->SetBinContent(i+1, yvalue); h1->SetBinError(i+1, yerr);
  }

  h1->Draw("");

  Tcan->Update();


  Tcan->cd();
  Tcan->SaveAs("/afs/cern.ch/work/p/pdas/www/emulator_phase2/12_3_0_pre4/singlepion/bugfix/response.pdf");
  Tcan->SaveAs("/afs/cern.ch/work/p/pdas/www/emulator_phase2/12_3_0_pre4/singlepion/bugfix/response.png");

  Tcan->Close();
  delete Tcan;


  TFile *MyFile = new TFile("response.root","NEW");
  h1->Write();
  MyFile->cd();
  MyFile->Write();
  delete h1; delete  MyFile; 
}
