for range in 5 10 15 20 25 30 35 40 45 50 55 60 70 80 90 100 120 140 160; do

cat>test1.C<<EOF
#include <iostream>
#include <fstream>
#include <set>
#include<TFile.h>
#include<TTree.h>
#include<TH1.h>
#include<TH2.h>
#include<THStack.h>
#include<TGraphErrors.h>
#include<TGraphAsymmErrors.h>
#include<TCanvas.h>
#include<TFrame.h>
#include<TLegend.h>
#include<vector>
#include<iostream>
#include<TMath.h>
#include<TROOT.h>
#include<TInterpreter.h>
#include<TStyle.h>
#include<TChain.h>
#include<TString.h>
#include<TPaveStats.h>
#include<TPad.h>
#include<TLatex.h>
#include "TEfficiency.h"
#include "TAxis.h"
#include <RooHistPdf.h>
#include <RooFormulaVar.h>
#include <RooPlot.h>
#include <RooAddPdf.h>
#include <RooFitResult.h>
#include <RooDataHist.h>
#include <RooDataSet.h>
#include <RooFit.h>
#include <RooArgList.h>
#include <RooRealVar.h>
#include <RooArgSet.h>
#include <RooLinkedList.h>
#include <RooVoigtian.h>
using namespace RooFit;


void test1(){

  TCanvas* c1 = new TCanvas("c1","", 100, 20, 1000, 800);
  c1->Range(0,0,1,1);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetFrameBorderMode(0);
  c1->Draw();

  gStyle->SetOptStat(0);
  TFile *fdata1 = TFile::Open("hfile_calib1.root");
  TH1F  *h1 = (TH1F*)fdata1->Get("reso${range}");

  double histomean = h1->GetMean();
  double histomean_err = h1->GetMeanError();

  //fstream out1;
  //out1.open("scale_pfcluster_doubleelectron_nofit.txt", fstream::in | fstream::out | fstream::app);
  //out1<<"reso${range}"<<"\t"<<histomean<<"\t"<<histomean_err<<endl;
  //out1.close();

  double xlow = h1->GetBinLowEdge(1);
  double xhigh = h1->GetBinLowEdge(h1->GetNbinsX())+h1->GetBinWidth(h1->GetNbinsX());
  RooRealVar x("x", "x", xlow, xhigh);
  RooDataHist data_hist("data_hist", "data_hist", RooArgList(x), h1);
  double m_da = h1->GetMean();
  double um_da = h1->GetMean() - h1->GetRMS();
  double uM_da = h1->GetMean() + h1->GetRMS();

  double mean = 1.;
  double error = 0.;
  RooRealVar v_m("v_m", "v_m", 0, 0., 2.);
  v_m.setVal(data_hist.mean(x));
  v_m.setRange(data_hist.mean(x) - data_hist.sigma(x), data_hist.mean(x) + data_hist.sigma(x));
  RooRealVar gamma_Z0("gamma_Z0_U", "Z0 width", 2.3, 0., 5., "GeV"); // gamma
  RooRealVar g_w("g_w", "width Gaus", 1., 0. , 2., "GeV"); // sigma
  RooVoigtian voigt("voigt", "Voigtian", x, v_m, gamma_Z0, g_w);
  RooPlot* frame = x.frame(Title("reso${range}"));
  frame->GetXaxis()->SetTitle("L1 p_{T} / Gen p_{T}");
  frame->GetYaxis()->SetTitleOffset(1.5);
  data_hist.plotOn(frame);

  RooFitResult *result = voigt.fitTo(data_hist, RooFit::Minimizer("Minuit","Migrad"), RooFit::Strategy(2), RooFit::SumW2Error(0), RooFit::Save(1), RooFit::PrintLevel(-1));
  voigt.plotOn(frame, RooFit::Components("voigt"), RooFit::VisualizeError(*result,1), RooFit::FillColor(kGray));
  voigt.plotOn(frame, RooFit::LineColor(kGray));
  mean = v_m.getVal();
  error = v_m.getError();

  if(frame->chiSquare() < 100.){
    frame->Draw();
  }
  else{
    data_hist.Draw();
  }

  fstream out2;
  out2.open("scale_pfcluster_calib1_fitted.txt", fstream::in | fstream::out | fstream::app);
  out2<<"reso${range}"<<"\t"<<mean<<"\t"<<error<<endl;
  out2.close();

  TLatex *latex = new TLatex();
  latex->SetNDC();
  latex->SetTextFont(42);
  latex->SetTextSize(0.04);
  latex->SetTextColor(kBlack);
  latex->DrawLatex(0.6, 0.840, Form("#scale[0.8]{Histogram mean = %.3f}", histomean));
  latex->DrawLatex(0.6, 0.780, Form("#scale[0.8]{Fitted mean = %.3f}", mean));
  latex->DrawLatex(0.6, 0.720, Form("#scale[0.8]{Fit #chi^{2} = %.3f}", frame->chiSquare()));

  c1->SaveAs("/afs/cern.ch/work/p/pdas/www/emulator_phase2/12_3_0_pre4/singlepion/calib/reso${range}_1.pdf");
  c1->SaveAs("/afs/cern.ch/work/p/pdas/www/emulator_phase2/12_3_0_pre4/singlepion/calib/reso${range}_1.png");

}

EOF

root -l -b -q test1.C

rm -f test1.C

done
