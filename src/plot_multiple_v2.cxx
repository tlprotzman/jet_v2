#include "TROOT.h"

#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1I.h"

#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"


int plot_multiple_v2() {
    TFile *r02_hc_file = new TFile("v2_hardcore_2.root");
    TFile *r04_hc_file = new TFile("v2_hardcore_4.root");
    TFile *r02_all_file = new TFile("v2_all_2.root");
    TFile *r04_all_file = new TFile("v2_all_4.root");

    TGraph *r02_hc_graph = r02_hc_file->Get<TGraphErrors>("Graph");
    TGraph *r04_hc_graph = r04_hc_file->Get<TGraphErrors>("Graph");
    TGraph *r02_all_graph = r02_all_file->Get<TGraphErrors>("Graph");
    TGraph *r04_all_graph = r04_all_file->Get<TGraphErrors>("Graph");

    TCanvas *c = new TCanvas("", "", 1000, 1000);
    c->SetMargin(0.15, 0.1, 0.15, 0.1);
    gStyle->SetOptStat(0);

    TH1I *dummy = new TH1I("", "", 1, 0, 40);
    dummy->GetYaxis()->SetRangeUser(-0.3, 0.4);
    dummy->GetXaxis()->SetTitle("Jet p_{T}");
    dummy->GetYaxis()->SetTitle("v_{2}^{ch}");
    // dummy->SetTitle("Jet v_{2}^{#hbox{ch}}");

    for (uint32_t i = 0; i < 2; i++) {
        r02_hc_graph->RemovePoint(0);
        r04_hc_graph->RemovePoint(0);
        r02_all_graph->RemovePoint(0);
        r04_all_graph->RemovePoint(0);
    }

    r02_hc_graph->SetMarkerColor(kRed);
    r04_hc_graph->SetMarkerColor(kBlue);
    r02_all_graph->SetMarkerColor(kRed);
    r04_all_graph->SetMarkerColor(kBlue);

    r02_hc_graph->SetLineColor(kRed);
    r04_hc_graph->SetLineColor(kBlue);
    r02_all_graph->SetLineColor(kRed+2);
    r04_all_graph->SetLineColor(kBlue+2);

    r02_hc_graph->SetMarkerStyle(kOpenSquare);
    r04_hc_graph->SetMarkerStyle(kOpenCircle);
    r02_all_graph->SetMarkerStyle(kFullSquare);
    r04_all_graph->SetMarkerStyle(kFullCircle);

    r02_hc_graph->SetMarkerSize(2);
    r04_hc_graph->SetMarkerSize(2);
    r02_all_graph->SetMarkerSize(2);
    r04_all_graph->SetMarkerSize(2);

    TLegend *legend = new TLegend(0.15, 0.15, 0.5, 0.35);
    legend->AddEntry(r02_hc_graph, "R=0.2, Hardcore Matched");
    legend->AddEntry(r02_all_graph, "R=0.2, Hardcore Only");
    legend->AddEntry(r04_hc_graph, "R=0.4, Hardcore Matched");
    legend->AddEntry(r04_all_graph, "R=0.4, Hardcore Only");

    dummy->Draw("AXIS");
    r02_hc_graph->Draw("p");
    r04_hc_graph->Draw("p");
    r02_all_graph->Draw("p");
    r04_all_graph->Draw("p");

    legend->Draw();

    TLatex *text = new TLatex();
    text->SetTextSize(0.03);
    text->DrawLatexNDC(0.18, 0.46, "STAR Internal");
    // text->DrawLatexNDC(0.18, 0.46, "Jet v_{2}^{#hbox{ch}}");
    text->DrawLatexNDC(0.18, 0.42, "#sqrt{S_{NN}}=200 GeV Ru+Ru/Zr+Zr");
    text->DrawLatexNDC(0.18, 0.38, "20-40\% Central");

    c->SaveAs("plots/v2_multiple.png");
    c->SaveAs("plots/v2_multiple.c");

    return 0;
}