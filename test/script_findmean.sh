for range in 5 10 15 20 25 30 35 40 45 50 55 60 70 80 90 100 120 140 160; do

cat>test1.C<<EOF
{
  TFile *fdata1 = TFile::Open("hfile_5x5.root");
  TH1F  *h1 = (TH1F*)fdata1->Get("reso${range}");

  double histomean = h1->GetMean();
  double histomean_err = h1->GetMeanError();

  fstream out;
  out.open("scale_pfcluster_5x5.txt", fstream::in | fstream::out | fstream::app);
  out<<"reso${range}"<<"\t"<<histomean<<"\t"<<histomean_err<<endl;
  out.close();
}
EOF

root -l -b -q test1.C

rm -f test1.C

done

cat>test2.C<<EOF
{
  TCanvas* c1 = new TCanvas("c1","", 100, 20, 1000, 800);
  c1->Range(0,0,1,1);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetFrameBorderMode(0);
  c1->Draw();


  TLine *line = new TLine(0, 1., 200, 1.);
  line->SetLineColor(kBlack);
  line->SetLineWidth(1.);
  line->SetLineStyle(3);


  //Tcan->SetGrid();

  //TLatex *latex = new TLatex();
  //latex->SetNDC();
  //latex->SetTextFont(42);
  //latex->SetTextColor(kBlack);

  //Tcan->cd();     /* Set current canvas */
  //Tcan->SetFillColor(0);

  ifstream f1;
  f1.open("scale_pfcluster.txt");
  string s;
  float yvalue, yerr;
  Float_t bins[] = {5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100, 120, 140, 160, 200};
  Int_t  binnum = 19;
  TH1F* h1 = new TH1F("h1","h1", binnum, bins);

  for(int i = 0; i < binnum; i++){
    f1>>s>>yvalue>>yerr;
    h1->SetBinContent(i+1, yvalue); h1->SetBinError(i+1, yerr);
  }


  TGraphErrors* gr1 = new TGraphErrors(h1);
  gr1->SetTitle("");
  gr1->SetMarkerStyle(8);
  gr1->SetMarkerSize(0.8);
  gr1->SetMarkerColor(kRed);
  gr1->SetLineColor(kRed);
  gr1->SetLineWidth(2);

  ifstream f2;
  f2.open("scale_pfcluster_5x5.txt");
  string s1;
  float yvalue1, yerr1;
  TH1F* h2 = new TH1F("h2","h2", binnum, bins);

  for(int i = 0; i < binnum; i++){
    f2>>s1>>yvalue1>>yerr1;
    h2->SetBinContent(i+1, yvalue1); h2->SetBinError(i+1, yerr1);
  }


  TGraphErrors* gr2 = new TGraphErrors(h2);
  gr2->SetTitle("");
  gr2->SetMarkerStyle(8);
  gr2->SetMarkerSize(0.8);
  gr2->SetMarkerColor(kBlue);
  gr2->SetLineColor(kBlue);
  gr2->SetLineWidth(2);

  gr1->GetXaxis()->SetTitle("Gen p_{T} [GeV]");
  gr1->GetXaxis()->SetRangeUser(0, 200);
  gr1->GetXaxis()->SetTitleSize(0.05);
  gr1->GetXaxis()->SetTitleOffset(0.9);

  gr1->GetYaxis()->SetTitle("<Reco p_{T}/Gen p_{T}>");
  gr1->GetYaxis()->SetRangeUser(0., 1.5);
  gr1->GetYaxis()->SetTitleSize(0.05);
  gr1->GetYaxis()->SetTitleOffset(0.9);

  gr1->Draw("apz");
  gr2->Draw("pez same");
  line->Draw("same");


  TLegend *legend1 = new TLegend(0.2, 0.7, 0.4, 0.85);
  legend1->SetTextFont(42);
  legend1->SetLineColor(0);
  legend1->SetTextSize(0.04);
  legend1->SetFillColor(0);
  legend1->AddEntry(gr1, "PF cluster", "lp");
  legend1->AddEntry(gr2, "5x5 cluster", "lp");
  legend1->Draw("same");

  TLatex *t2a = new TLatex(0.5,0.91," #bf{CMS} #it{Phase-2 GCT emulator}                         0 PU ");
  t2a->SetNDC();
  t2a->SetTextFont(42);
  t2a->SetTextSize(0.05);
  t2a->SetTextAlign(20);
  t2a->Draw("same");

  TLatex *latex = new TLatex();
  latex->SetNDC();
  latex->SetTextFont(42);
  latex->SetTextColor(kBlack);
  latex->DrawLatex(0.54, 0.840, "#scale[0.8]{EG Barrel}");
  latex->DrawLatex(0.54, 0.780, "#scale[0.8]{SinglePion Pt 0 to 200}");
  latex->DrawLatex(0.54, 0.720, "#scale[0.8]{Gen p_{T} > 5, |#eta^{Gen}| < 1.4841}");

  //c1->Update();

  //Tcan->cd();
  c1->SaveAs("/afs/cern.ch/work/p/pdas/www/emulator_phase2/12_3_0_pre4/singlepion/bugfix/response_5x5.pdf");
  c1->SaveAs("/afs/cern.ch/work/p/pdas/www/emulator_phase2/12_3_0_pre4/singlepion/bugfix/response_5x5.png");

  c1->Close();
  delete c1;

  TFile *MyFile = new TFile("response_5x5.root","RECREATE");
  h2->Write();
  MyFile->cd();
  MyFile->Write();
  delete h1; delete gr1; delete  MyFile;
  delete h2; delete gr2;
}
EOF

root -l -b -q test2.C

rm -rf test2.C

#rm -rf scale_pfcluster.txt

