#include "../baseCodeForPlots/tdrstyle.C"
#include "../baseCodeForPlots/comparisonPlots.C"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TH1F.h"

#include "TAxis.h"
#include "TChain.h"

#include "THStack.h"
#include "TLegend.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TTree.h"

#include <assert.h>
#include <cmath>
#include <iostream>
#include <math.h>
#include <string>
#include <TMath.h>
#include <vector>
void plot_poster(){
  setTDRStyle();
  TCanvas* c1 = new TCanvas("c1","", 100, 20, 1000, 800);
  c1->Range(0,0,1,1);
  c1->SetFillColor(0);
  //c1->SetBorderMode(0);
  //c1->SetBorderSize(2);
  c1->SetFrameBorderMode(0);
  c1->Draw();

  TLine *line = new TLine(0, 1., 200, 1.);
  line->SetLineColor(kBlack);
  line->SetLineWidth(1.);
  line->SetLineStyle(3);

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

  gr1->GetXaxis()->SetTitle("Gen p_{T} [GeV]");
  gr1->GetXaxis()->SetRangeUser(0, 200);
  gr1->GetXaxis()->SetTitleSize(0.06);
  //gr1->GetXaxis()->SetTitleOffset(0.6);

  gr1->GetYaxis()->SetTitle("<L1 p_{T}/Gen p_{T}>");
  gr1->GetYaxis()->SetRangeUser(0., 1.5);
  gr1->GetYaxis()->SetTitleSize(0.06);
  //gr1->GetYaxis()->SetTitleOffset(0.6);

  gr1->Draw("apz");
  line->Draw("same");

  TLegend *legend1 = new TLegend(0.2, 0.7, 0.4, 0.85);
  legend1->SetTextFont(42);
  legend1->SetLineColor(0);
  legend1->SetTextSize(0.04);
  legend1->SetFillColor(0);
  legend1->AddEntry(gr1, "PF cluster", "lp");
  legend1->Draw("same");

//  TLatex *t2a = new TLatex(0.5,0.96," #bf{CMS} #it{Phase-2 GCT simulation}                       14 TeV (0 PU) ");
//  t2a->SetNDC();
//  t2a->SetTextFont(42);
//  t2a->SetTextSize(0.04);
//  t2a->SetTextAlign(20);
//  t2a->Draw("same");

  TLatex *latex = new TLatex();
  latex->SetNDC();
  latex->SetTextFont(42);
  latex->SetTextColor(kBlack);
  latex->DrawLatex(0.16, 0.960, "#scale[1.0]{#bf{CMS}} #scale[0.8]{#it{Phase 2 GCT simulation}}");
  latex->DrawLatex(0.79, 0.960, "#scale[0.8]{14 TeV (0 PU)}");
  latex->DrawLatex(0.2, 0.840, "#scale[0.8]{Work in progress}");
  latex->DrawLatex(0.62, 0.900, "#scale[0.8]{EG Barrel}");
  latex->DrawLatex(0.62, 0.840, "#scale[0.8]{SinglePion Pt 0 to 200}");
  latex->DrawLatex(0.62, 0.780, "#scale[0.8]{Gen p_{T} > 5, |#eta^{Gen}| < 1.4841}");

  c1->SaveAs("/afs/cern.ch/work/p/pdas/www/emulator_phase2/12_3_0_pre4/singlepion/response.pdf");
  c1->SaveAs("/afs/cern.ch/work/p/pdas/www/emulator_phase2/12_3_0_pre4/singlepion/response.png");

  c1->Close();
  delete c1;

  TFile *MyFile = new TFile("response.root","RECREATE");
  h1->Write();
  MyFile->cd();
  MyFile->Write();
  delete h1; delete gr1; delete  MyFile;
}
