for range in 5 10 15 20 25 30 35 40 45 50 55 60 70 80 90 100 120 140 160; do

cat>test1.C<<EOF
{

  TCanvas* c1 = new TCanvas("c1","", 100, 20, 1000, 800);
  c1->Range(0,0,1,1);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetFrameBorderMode(0);
  c1->Draw();

  gStyle->SetOptStat(0);
  TFile *fdata1 = TFile::Open("hfile.root");
  TH1F  *h1 = (TH1F*)fdata1->Get("reso${range}");

  h1->Draw("hist");

  double histomean = h1->GetMean();
  double histomean_err = h1->GetMeanError();

  std::stringstream ss, ss1;
  ss << histomean; 
  ss1 << histomean_err;
  const char* str = ss.str().c_str();
  const char* str1 = ss1.str().c_str();

  TLatex *test1 = new TLatex(0.6, 0.8, "Mean =");
  test1->SetNDC();
  test1->Draw("same");
  TLatex *test2 = new TLatex(0.68, 0.8, str);
  test2->SetNDC();
  test2->Draw("same");

  TLatex *test3 = new TLatex(0.6, 0.7, "Error =");
  test3->SetNDC();
  test3->Draw("same");
  TLatex *test4 = new TLatex(0.68, 0.7, str1);
  test4->SetNDC();
  test4->Draw("same");

  c1->SaveAs("/afs/cern.ch/work/p/pdas/www/emulator_phase2/12_3_0_pre4/singlepion/bugfix/reso${range}.pdf");
  c1->SaveAs("/afs/cern.ch/work/p/pdas/www/emulator_phase2/12_3_0_pre4/singlepion/bugfix/reso${range}.png");

}
EOF

root -l -b -q test1.C

rm -f test1.C

done

