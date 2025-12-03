#include "MHist.h"
#include "MRootGraphic.h"
#include "MRootIO.h"
#include "TGraph.h"
#include "TLegend.h"

TFile *file_input = nullptr;
TFile *file_output = nullptr;

void DrawHist(TString name_observable = "", TString name_ncls = "",
              TString particle = "") {
  TString name_hist =
      "fEta_fPhi_" + name_observable + "_" + name_ncls + "_" + particle;

  auto hist = (THnD *)file_input->Get(name_hist);
  MHnTool hnTool(hist);

  vector<TH1D *> vec_eta_diff;
  for (int i = 0; i < hnTool.GetNbins(0); i++) {
    double eta_low = hnTool.hN->GetAxis(0)->GetBinLowEdge(i + 1);
    double eta_high = eta_low + hnTool.hN->GetAxis(0)->GetBinWidth(i + 1);
    TH2D *obs_eta_diff = hnTool.Project(2, 1, {i + 1});
    TString cond(hist->GetTitle());
    TString cond_etaRange = Form("%.2f < #eta < %.2f", eta_low, eta_high);
    cond += ", " + cond_etaRange;
    obs_eta_diff->SetTitle(cond);
    MRootGraphic::StyleHistCommon(obs_eta_diff);
    obs_eta_diff->SetDirectory(file_output);
    TH1D *mean_obs = obs_eta_diff->ProfileX(
        Form("mean_%s_eta_%d", name_observable.Data(), i));
    mean_obs->GetYaxis()->SetTitle("#LT " + name_observable + " #GT");
    mean_obs->SetTitle(cond + ", Profile of " + name_observable);
    mean_obs->SetDirectory(file_output);
    MRootGraphic::StyleHistCommon(mean_obs);
  }
}

void PhiCheck(
    TString path_input = "/home/szhu/work/alice/analysis/pid_qa/test/test.root",
    TString path_output =
        "/home/szhu/work/alice/analysis/pid_qa/test/output_test.root") {
  // TFile *file = new TFile(path_input, "READ");
  // TFile *file_output = new TFile(path_output, "RECREATE");
  file_input = new TFile(path_input, "READ");
  file_output = new TFile(path_output, "RECREATE");
  for (auto particle : {"Electron", "Pion", "Kaon", "Proton"}) {
    for (auto ncls : {"highNcls", "lowNcls", "AllNcls"}) {
      for (auto observable : {"dEdx", "fNSigTPC"}) {
        DrawHist(observable, ncls, particle);
      }
    }
  }
  file_output->Write();
  file_output->Close();
}