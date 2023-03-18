#include "TROOT.h"

#include "TFile.h"
#include "TDirectory.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TH1I.h"

#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TBox.h"

#include <iostream>
#include <vector>

const int text_font = 42;
const double text_size = 0.03;

bool large = false;
bool small = false;
bool systems = false;
bool trackeff = true;


int plot_multiple_v2(std::string which="") {
    if (which != "") {
        large = which == "large";
        small = which == "small";
        systems = which == "systems";
        trackeff = which == "trackeff";
    }
    TFile *r02_hc_file = new TFile("testout.root");
    TFile *r02_old_file = new TFile("v2_hardcore_2.root");
    TFile *r02_tracking_eff_file = new TFile("testout.trackingefficiency.root");
    // TFile *r02_all_file = new TFile("testout.root");
    // TFile *r04_all_file = new TFile("v2_all_4.root");

    TFile *alice_result = new TFile("alice_v2_result.root");

    TGraph *r02_hc_graph = r02_hc_file->Get<TGraphErrors>("v2_hc");
    TGraph *r02_avg_ep = r02_hc_file->Get<TGraphErrors>("v2_hc");
    TGraph *r02_all_graph = r02_hc_file->Get<TGraphErrors>("v2_all");
    TGraph *r02_no_mod = r02_hc_file->Get<TGraphErrors>("v2_no_mod");
    TGraph *r02_extra_mod = r02_hc_file->Get<TGraphErrors>("v2_extra_mod");
    TGraph *r02_hc_only = r02_hc_file->Get<TGraphErrors>("v2_hc_only");
    TGraph *r02_trackeff = r02_tracking_eff_file->Get<TGraphErrors>("v2_hc");
    // TGraph *r04_all_graph = r04_all_file->Get<TGraphErrors>("Graph");
    TDirectory *alice_dir = alice_result->GetDirectory("Table 2");
    TH1 *alice_v2_graph = alice_dir->Get<TH1F>("Hist1D_y1");
    TH1 *alice_v2_stat_error = alice_dir->Get<TH1F>("Hist1D_y1_e1");
    TH1 *alice_v2_shape_error = alice_dir->Get<TH1F>("Hist1D_y1_e2");
    TH1 *alice_v2_corr_error_plus = alice_dir->Get<TH1F>("Hist1D_y1_e3plus");
    TH1 *alice_v2_corr_error_minus = alice_dir->Get<TH1F>("Hist1D_y1_e3minus");

    TGraphErrors *r02_systematics = new TGraphErrors();
    std::vector<double> point_systematics;
    for (uint i = 0; i < r02_hc_graph->GetN(); i++) {
        std::cerr << i << std::endl;
        double point_x, point_a_x, point_b_x;
        double point_y, point_a_y, point_b_y;
        double ex, ey;
        r02_hc_graph->GetPoint(i, point_x, point_y);
        r02_no_mod->GetPoint(i, point_a_x, point_a_y);
        r02_extra_mod->GetPoint(i, point_b_x, point_b_y);
        double average_delta = 0.5 * abs((point_a_y - point_y) + abs(point_b_y - point_y));
        r02_systematics->SetPoint(i, point_x, point_y);
        r02_systematics->SetPointError(i, r02_hc_graph->GetErrorX(i), average_delta);
        point_systematics.push_back(average_delta);
    }
    // Manually tweak the 10-12.5 and 15-17.5 bin systematics
    r02_systematics->SetPointError(2, r02_hc_graph->GetErrorX(2), point_systematics[3]);
    r02_systematics->SetPointError(4, r02_hc_graph->GetErrorX(4), 0.5 * (point_systematics[3] + point_systematics[5]));

    TH1 *alice_v2_shape_error_plot = (TH1*)alice_v2_graph->Clone("alice_shape_err");
    TGraphAsymmErrors *alice_v2_corr_error_plot = new TGraphAsymmErrors(alice_v2_graph);
    for (uint i = 1; i <= alice_v2_graph->GetNbinsX(); i++) {
        alice_v2_graph->SetBinError(i, alice_v2_stat_error->GetBinContent(i));
        alice_v2_shape_error_plot->SetBinError(i, alice_v2_shape_error->GetBinContent(i));
        alice_v2_shape_error_plot->SetBinError(i, alice_v2_shape_error->GetBinContent(i));

        alice_v2_corr_error_plot->SetPoint(i-1, alice_v2_graph->GetBinCenter(i), alice_v2_graph->GetBinContent(i));
        alice_v2_corr_error_plot->SetPointEXlow(i-1, alice_v2_graph->GetBinWidth(i) / 2.);
        alice_v2_corr_error_plot->SetPointEXhigh(i-1, alice_v2_graph->GetBinWidth(i) / 2.);
        std::cout << "up: " << alice_v2_corr_error_plus->GetBinContent(i) << " low: " << alice_v2_corr_error_minus->GetBinContent(i) << std::endl;
        alice_v2_corr_error_plot->SetPointEYhigh(i-1, alice_v2_corr_error_plus->GetBinContent(i));
        alice_v2_corr_error_plot->SetPointEYlow(i-1, abs(alice_v2_corr_error_minus->GetBinContent(i)));
    }

    alice_v2_corr_error_plot->SetMarkerSize(0.);
    alice_v2_corr_error_plot->SetFillStyle(1001);
    alice_v2_corr_error_plot->SetMarkerColor(kGray);
    alice_v2_corr_error_plot->SetFillColor(18);
    alice_v2_corr_error_plot->SetLineWidth(0);

    alice_v2_shape_error_plot->SetMarkerSize(0.);
    alice_v2_shape_error_plot->SetLineWidth(2);
    alice_v2_shape_error_plot->SetLineColor(kCyan-8);
    alice_v2_shape_error_plot->SetFillStyle(0);

    r02_systematics->SetMarkerSize(0);
    r02_systematics->SetFillStyle(0);
    r02_systematics->SetFillColor(kGray);
    r02_systematics->SetLineWidth(2);
    r02_systematics->SetLineColor(kBlack);




    TCanvas *c = new TCanvas("", "", 1400, 1000);
    c->SetMargin(0.1, 0.02, 0.1, 0.05);
    gStyle->SetOptStat(0);

    TH1I *dummy;
    TBox *alice_shape_error_box;
    TBox *r02_systematics_box;
    if (large) {
        dummy = new TH1I("", "", 1, 5, 92);
        alice_shape_error_box = new TBox(38.3, 0.1655, 45.65, 0.1755); // for 5 to 92
        r02_systematics_box = new TBox(38.3, 0.191, 45.65, 0.201); // for 5 to 92
        dummy->GetYaxis()->SetRangeUser(-0.01, 0.22);
    }
    else if (small) {
        dummy = new TH1I("", "", 1, 5, 42);
        alice_shape_error_box = new TBox(19.15, 0.1655, 22.27, 0.1755); // for 5 to 42
        r02_systematics_box = new TBox(19.15, 0.191, 22.27, 0.201); // for 5 to 42
        dummy->GetYaxis()->SetRangeUser(-0.01, 0.22);
    }
    else if (systems || trackeff) {
        dummy = new TH1I("", "", 1, 7, 24);
        alice_shape_error_box = new TBox(38.3, 0.17, 45.65, 0.18); // Not used
        r02_systematics_box = new TBox(38.3, 0.17, 45.65, 0.18); // Not used
        dummy->GetYaxis()->SetRangeUser(-0.2, 0.35);
    }
    dummy->GetXaxis()->SetTitle("p_{T, jet}^{reco, ch}, p_{T, jet}^{ch} (GeV/c)");
    dummy->GetXaxis()->SetTitleSize(0.042);
    dummy->GetXaxis()->SetTitleOffset(1.02);
    dummy->GetYaxis()->SetTitle("v_{2}^{ch jet}");
    dummy->GetYaxis()->SetTitleSize(0.05);
    dummy->GetYaxis()->SetTitleOffset(0.9);
    // dummy->SetTitle("Jet v_{2}^{#hbox{ch}}");

    alice_shape_error_box->SetFillColor(kWhite);
    alice_shape_error_box->SetLineColor(kCyan+2);
    alice_shape_error_box->SetLineWidth(2);

    r02_systematics_box->SetFillColor(kWhite);
    r02_systematics_box->SetLineColor(kBlack);
    r02_systematics_box->SetLineWidth(2);

    for (uint32_t i = 0; i < r02_all_graph->GetN(); i++) {
        r02_all_graph->SetPointX(i, r02_all_graph->GetPointX(i) - 0.2);
        r02_hc_only->SetPointX(i, r02_hc_graph->GetPointX(i) + 0.2);
        r02_no_mod->SetPointX(i, r02_hc_graph->GetPointX(i) - 0.2);
        r02_extra_mod->SetPointX(i, r02_hc_graph->GetPointX(i) + 0.2);
        r02_trackeff->SetPointX(i, r02_hc_graph->GetPointX(i) + 0.2);
    }

    // for (uint32_t i = 0; i < 2; i++) {
        r02_hc_graph->RemovePoint(0);
        r02_hc_graph->RemovePoint(0);
        r02_extra_mod->RemovePoint(0);
        r02_extra_mod->RemovePoint(0);
        r02_no_mod->RemovePoint(0);
        r02_no_mod->RemovePoint(0);
        r02_trackeff->RemovePoint(0);
        r02_trackeff->RemovePoint(0);
        // r02_avg_ep->RemovePoint(0);
        r02_all_graph->RemovePoint(0);
        r02_all_graph->RemovePoint(0);
        r02_hc_only->RemovePoint(0);
        r02_hc_only->RemovePoint(0);
        r02_systematics->RemovePoint(0);
        r02_systematics->RemovePoint(0);
    // }

    for (uint32_t i = 0; i < 3; i++) {
        r02_hc_graph->RemovePoint(5);
        r02_no_mod->RemovePoint(5);
        r02_extra_mod->RemovePoint(5);
        r02_avg_ep->RemovePoint(5);
        r02_all_graph->RemovePoint(5);
        r02_hc_only->RemovePoint(5);
        r02_systematics->RemovePoint(5);
        r02_trackeff->RemovePoint(5);

    }

    r02_hc_graph->SetMarkerColor(kRed);
    r02_no_mod->SetMarkerColor(kViolet);
    r02_extra_mod->SetMarkerColor(kTeal+2);
    r02_avg_ep->SetMarkerColor(kBlue);
    r02_all_graph->SetMarkerColor(kGreen+2);
    r02_hc_only->SetMarkerColor(kBlue);
    alice_v2_graph->SetMarkerColor(kGreen-5);
    r02_trackeff->SetMarkerColor(kOrange+1);

    r02_hc_graph->SetLineColor(kRed);
    r02_no_mod->SetLineColor(kMagenta);
    r02_extra_mod->SetLineColor(kTeal+2);
    r02_avg_ep->SetLineColor(kBlue);
    r02_all_graph->SetLineColor(kGreen+4);
    r02_hc_only->SetLineColor(kBlue+2);
    alice_v2_graph->SetLineColor(kGreen-5);
    r02_trackeff->SetLineColor(kOrange+1);

    r02_hc_graph->SetMarkerStyle(kFullSquare);
    r02_no_mod->SetMarkerStyle(kFullDiamond);
    r02_extra_mod->SetMarkerStyle(kFullCrossX);
    r02_avg_ep->SetMarkerStyle(kOpenCircle);
    r02_all_graph->SetMarkerStyle(kFullCrossX);
    r02_hc_only->SetMarkerStyle(kFullCircle);
    alice_v2_graph->SetMarkerStyle(kOpenCircle);
    r02_trackeff->SetMarkerStyle(kFullCircle);

    r02_hc_graph->SetMarkerSize(2);
    r02_no_mod->SetMarkerSize(2);
    r02_extra_mod->SetMarkerSize(2);
    r02_avg_ep->SetMarkerSize(2);
    r02_all_graph->SetMarkerSize(2);
    r02_hc_only->SetMarkerSize(2);
    alice_v2_graph->SetMarkerSize(2);
    r02_trackeff->SetMarkerSize(2);

    r02_hc_graph->SetLineWidth(2);
    r02_no_mod->SetLineWidth(2);
    r02_extra_mod->SetLineWidth(2);
    r02_avg_ep->SetLineWidth(2);
    r02_all_graph->SetLineWidth(2);
    alice_v2_graph->SetLineWidth(2);
    r02_trackeff->SetLineWidth(2);

    TLegend *legend;
    if (large || small) {
        legend = new TLegend(0.42, 0.69, 0.85, 0.94);
        legend->AddEntry(r02_hc_graph, "R=0.2 Matched, HC Constituents > 2 GeV");
        legend->AddEntry(r02_systematics, "Background Modulation Systematic Uncertainty");
        legend->AddEntry(alice_v2_graph, "ALICE R=0.2 Pb+Pb #sqrt{s_{NN}}=2.76 TeV, 30-50\%");
        legend->AddEntry(alice_v2_shape_error_plot, "ALICE Systematic Uncertainty (Shape)");
        legend->AddEntry(alice_v2_corr_error_plot, "ALICE Systematic Uncertainty (Correlated)");
    }
    else if (systems) {
        legend = new TLegend(0.42, 0.71, 0.85, 0.91);
        legend->AddEntry(r02_all_graph, "R=0.2 All Jets, Statistical Error Only");
        legend->AddEntry(r02_hc_only, "R=0.2 Hard Cores Only, Statistical Error Only");
        legend->AddEntry(r02_hc_graph, "R=0.2 Matched, HC Constituents > 2 GeV");
    }
    else if (trackeff) {
        legend = new TLegend(0.42, 0.69, 0.85, 0.94);
        legend->AddEntry(r02_hc_graph, "R=0.2 Matched, HC Constituents > 2 GeV");
        legend->AddEntry(r02_systematics, "Background Modulation Systematic Uncertainty");
        legend->AddEntry(r02_trackeff, "R=0.2 Matched, 96% Track Efficiency");
    }
    legend->SetBorderSize(0);
    legend->SetTextSize(text_size);
    legend->SetTextFont(text_font);
    // legend->AddEntry(r02_avg_ep, "R=0.2, Average EP");
    // legend->AddEntry(r02_no_mod, "R=0.2, 2\% Charged Particle v_{2}");
    // legend->AddEntry(r02_extra_mod, "R=0.2, 6\% Charged Particle v_{2}");

    TLatex *text = new TLatex();
    text->SetTextSize(text_size);
    text->SetTextFont(text_font);

    dummy->Draw("AXIS");
    legend->Draw();
    if (large || small) {
        alice_v2_corr_error_plot->Draw("E2 same");
        alice_v2_shape_error_plot->Draw("E2 same");
        alice_v2_graph->Draw("e0 same");
        r02_systematics->Draw("E2 same");
        r02_hc_graph->Draw("p");

        alice_shape_error_box->Draw("l");
        r02_systematics_box->Draw("l");
        text->SetTextFont(52);
        text->DrawLatexNDC(0.526, 0.66, "Phys. Lett B 753 (2016) 511");
        text->SetTextFont(42);
    }
    else if (systems) {
        r02_systematics->Draw("E2 same");
        r02_hc_graph->Draw("p");
        r02_all_graph->Draw("p");
        r02_hc_only->Draw("p");
    }
    else if (trackeff) {
        r02_systematics->Draw("E2 same");
        r02_hc_graph->Draw("p");
        r02_trackeff->Draw("p");
    }
    // r02_no_mod->Draw("p");
    // r02_extra_mod->Draw("p");
    // r02_avg_ep->Draw("p");


    text->DrawLatexNDC(0.135, 0.91, "STAR Preliminary");
    // text->DrawLatexNDC(0.15, 0.46, "Jet v_{2}^{#hbox{ch}}");
    text->DrawLatexNDC(0.135, 0.86, "#sqrt{S_{NN}}=200 GeV Ru+Ru & Zr+Zr");
    text->DrawLatexNDC(0.135, 0.81, "Jet |#eta| < 0.8");
    text->DrawLatexNDC(0.135, 0.77, "Event Plane 2.1 < |#eta| < 5.1");
    text->DrawLatexNDC(0.135, 0.73, "20-60\% Mid-Central");
    text->DrawLatexNDC(0.135, 0.69, "No Correction on Jet p_{T} Applied");
    // text->DrawLatexNDC(0.135, 0.69, "Statistical Uncertainties Only");

    // c->SetLogx();

    if (large) {
        c->SaveAs("plots/v2_multiple_large.png");
        c->SaveAs("plots/v2_multiple_large.pdf");
        c->SaveAs("plots/v2_multiple_large.c");
    }
    if (small) {
        c->SaveAs("plots/v2_multiple_small.png");
        c->SaveAs("plots/v2_multiple_small.pdf");
        c->SaveAs("plots/v2_multiple_small.c");
    }
    if (systems) {
        c->SaveAs("plots/v2_multiple_systems.png");
        c->SaveAs("plots/v2_multiple_systems.pdf");
        c->SaveAs("plots/v2_multiple_systems.c");
    }
    if (trackeff) {
        c->SaveAs("plots/v2_multiple_trackeff.png");
        c->SaveAs("plots/v2_multiple_trackeff.pdf");
        c->SaveAs("plots/v2_multiple_trackeff.c");
    }

    return 0;
}