#include "MRootGraphic.h"
#include "MRootIO.h"
#include "TGraph.h"
#include "TLegend.h"

TGraph *GetGraph(TH1D *h1, TString name) {
  TGraph *g = new TGraph(h1->GetNbinsX());
  for (int i = 0; i < h1->GetNbinsX(); i++) {
    g->SetPoint(i, h1->GetBinCenter(i + 1), h1->GetBinContent(i + 1));
  }
  g->SetName(name);
  g->SetTitle(h1->GetTitle());
  g->GetXaxis()->SetTitle(h1->GetXaxis()->GetTitle());
  g->GetYaxis()->SetTitle(h1->GetYaxis()->GetTitle());
  return g;
}

void PlottingDetailed(TFile *file_nn, TFile *file_bb, TString path_output_graph,
                      TString tag_x, TString tag_ncls) {
  TString name_tag = tag_x + "_" + tag_ncls;
  TCanvas *c = new TCanvas("canvas_detailed_" + name_tag,
                           "canvas_detailed_" + name_tag, 800, 800);

#define GetGraphPlotlottingDetailed(nn_bb, mean_sigma, var)                    \
  TH1D *mean_sigma##_##nn_bb##_##var##elec = (TH1D *)file_##nn_bb->Get(        \
      tag_x + "_" + #var + "_" + tag_ncls + "_Electron_" #mean_sigma);         \
  TGraph *g_##mean_sigma##_##nn_bb##_##var##elec = GetGraph(                   \
      mean_sigma##_##nn_bb##_##var##elec,                                      \
      "g_" + tag_x + "_" + #var + "_" + tag_ncls + "_" + #nn_bb + "_elec");    \
  TH1D *mean_sigma##_##nn_bb##_##var##pion = (TH1D *)file_##nn_bb->Get(        \
      tag_x + "_" + #var + "_" + tag_ncls + "_Pion_" #mean_sigma);             \
  TGraph *g_##mean_sigma##_##nn_bb##_##var##pion = GetGraph(                   \
      mean_sigma##_##nn_bb##_##var##pion,                                      \
      "g_" + tag_x + "_" + #var + "_" + tag_ncls + "_" + #nn_bb + "_pion");

  GetGraphPlotlottingDetailed(nn, mean, dEdx);
  GetGraphPlotlottingDetailed(nn, sigma, delta_dEdx);
  GetGraphPlotlottingDetailed(bb, mean, dEdx);
  GetGraphPlotlottingDetailed(bb, sigma, delta_dEdx);

  TGraph *g_nn, *g_bb;
  c->Divide(2, 2);
#define SingleGraph(nn_bb, mean_sigma, var, species)                           \
  g_##mean_sigma##_##nn_bb##_##var##species

#define GroupGraph(...)                                                        \
  g_nn = SingleGraph(nn, __VA_ARGS__);                                         \
  g_bb = SingleGraph(bb, __VA_ARGS__);                                         \
  MRootGraphic::StyleHistCommonGraph(g_nn);                                    \
  MRootGraphic::StyleHistCommonGraph(g_bb);                                    \
  g_nn->SetLineColor(kRed);                                                    \
  g_bb->SetLineColor(kBlue);                                                   \
  g_nn->Draw("AL");                                                            \
  g_bb->Draw("L same");

  c->cd(1);
  GroupGraph(mean, dEdx, elec);
  c->cd(2);
  GroupGraph(sigma, delta_dEdx, elec);
  c->cd(3);
  GroupGraph(mean, dEdx, pion);
  c->cd(4);
  GroupGraph(sigma, delta_dEdx, pion);

  TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
  leg->SetFillStyle(0);
  leg->SetLineStyle(0);
  leg->SetTextSize(0.04);
  leg->SetLineColor(0);
  leg->AddEntry(g_nn, "NN", "L");
  leg->AddEntry(g_bb, "BB", "L");
  leg->Draw("same");

  c->SaveAs(path_output_graph + "_detailed_" + name_tag + ".json");
  c->Close();
  delete c;
}

void PlottingDetailed(TFile *file_nn, TFile *file_bb,
                      TString path_output_graph) {
  vector<TString> vec_tag_x = {"fTgl", "fFt0Occ"};
  vector<TString> vec_tag_ncls = {"highNcls", "lowNcls", "AllNcls"};

  for (const auto &tag_x : vec_tag_x)
    for (const auto &tag_ncls : vec_tag_ncls)
      PlottingDetailed(file_nn, file_bb, path_output_graph, tag_x, tag_ncls);
}

void PlottingSepPower(TFile *file_nn, TFile *file_bb,
                      TString path_output_graph) {
  vector<TString> vec_tag_x = {"fTgl", "fFt0Occ"};
  vector<TString> vec_tag_ncls = {"highNcls", "lowNcls", "AllNcls"};

  TCanvas *c_sepPower =
      new TCanvas("canvas_sepPower", "canvas_sepPower", 800, 900);
  c_sepPower->Divide(2, 3);

  int index = 1;
  for (const auto &tag_x : vec_tag_x) {
    for (const auto &tag_ncls : vec_tag_ncls) {
      c_sepPower->cd(index++);
      gPad->SetTopMargin(0.1);
      TString name_tag = tag_x + "_sepPower_" + tag_ncls;
      TH1D *sepPower_nn = (TH1D *)file_nn->Get(name_tag);
      TH1D *sepPower_bb = (TH1D *)file_bb->Get(name_tag);
      TGraph *g_sepPower_nn = GetGraph(sepPower_nn, "g_" + name_tag + "_nn");
      TGraph *g_sepPower_bb = GetGraph(sepPower_bb, "g_" + name_tag + "_bb");
      MRootGraphic::StyleHistCommonGraph(g_sepPower_nn);
      MRootGraphic::StyleHistCommonGraph(g_sepPower_bb);
      g_sepPower_nn->SetLineColor(kRed);
      g_sepPower_bb->SetLineColor(kBlue);
      double max =
          TMath::Max(sepPower_nn->GetMaximum(), sepPower_bb->GetMaximum());
      double min =
          TMath::Min(sepPower_nn->GetMinimum(), sepPower_bb->GetMinimum());
      double max_user = max + 0.1 * (max - min);
      double min_user = min - 0.1 * (max - min);
      g_sepPower_nn->GetYaxis()->SetRangeUser(min_user, max_user);
      g_sepPower_nn->Draw("AL");
      g_sepPower_bb->Draw("L same");

      if (index == 2) {
        TLegend *leg = new TLegend(0.5, 0.6, 0.7, 0.8);
        leg->SetFillStyle(0);
        leg->SetLineStyle(0);
        leg->SetLineColor(0);
        leg->SetTextSize(0.04);
        leg->AddEntry(g_sepPower_nn, "NN", "L");
        leg->AddEntry(g_sepPower_bb, "BB", "L");
        leg->Draw("same");
      }
    }
  }
  c_sepPower->SaveAs(path_output_graph + "_sepPower.json");
}

void SeparationPower(TString path_input_nn = "~/test/sepPower_nn.root",
                     TString path_input_bb = "~/test/sepPower_bb.root",
                     TString path_output_graph = "~/test/sepPower_graph") {
  TFile *file_nn = new TFile(path_input_nn, "READ");
  TFile *file_bb = new TFile(path_input_bb, "READ");

  MRootGraphic::StyleCommon();
  PlottingDetailed(file_nn, file_bb, path_output_graph);
  PlottingSepPower(file_nn, file_bb, path_output_graph);
}
