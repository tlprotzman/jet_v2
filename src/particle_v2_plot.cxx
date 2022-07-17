#include <TROOT.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>

#include <iostream>
#include <vector>

int colors[16] = {kAzure + 10, kAzure,
                         kViolet + 10, kViolet-3, kViolet - 9,
                         kPink + 10, kPink,
                         kOrange + 10, kOrange, kOrange - 9,
                         kSpring + 10, kSpring,
                         kTeal + 10, kTeal, kTeal - 9, kCyan + 4
                         };

char *cent_bins[8] = {"70-80%", "60-70%", "50-60%", "40-50%", "30-40%", "20-30%", "10-20%", "0-10%",};
int pt_bin_width  = 2;

void particle_v2_plot() {
    TFile infile("out.root");
    TH3 *v2_hist = infile.Get<TH3>("particle_v2");
    if (v2_hist == nullptr) {
        std::cout << "Could not find particle_v2" << std::endl;
        return;
    }

    std::vector<TGraph*> *v2_data = new std::vector<TGraph*>;

    TF1 *fourier_expansion = new TF1("fourier_expansion", "[0] * (1 + 2 * [1] * cos(2 * x))");
    fourier_expansion->SetParName(0, "scale");
    fourier_expansion->SetParName(1, "v2");

    std::vector<double> centrality;
    std::vector<double> centrality_err;
    std::vector<double> v2;
    std::vector<double> v2_err;

    for (int pt_bin = 1; pt_bin <= v2_hist->GetYaxis()->GetNbins(); pt_bin += pt_bin_width) {
        centrality.clear();
        v2.clear();
        v2_hist->GetYaxis()->SetRange(pt_bin, pt_bin);
        for (int cent_bin = 1; cent_bin <= v2_hist->GetZaxis()->GetNbins(); cent_bin += 2) {
            v2_hist->GetZaxis()->SetRange(cent_bin, cent_bin + 1);
            TH1 *dNdPhi = v2_hist->Project3D("x");
            dNdPhi->Fit(fourier_expansion, "q");
            centrality.push_back(v2_hist->GetZaxis()->GetBinCenter(cent_bin));
            std::cout << v2_hist->GetZaxis()->GetBinCenter(cent_bin) << std::endl;
            centrality_err.push_back(v2_hist->GetZaxis()->GetBinWidth(cent_bin));
            v2.push_back(fourier_expansion->GetParameter("v2"));
            v2_err.push_back(fourier_expansion->GetParError(fourier_expansion->GetParNumber("v2")));
        }
        TGraph *v2_graph = new TGraphErrors(centrality.size(), &centrality[0], &v2[0], &centrality_err[0], &v2_err[0]);
        v2_graph->SetTitle(Form("Centrality_Bin_%2i", pt_bin));
        v2_data->push_back(v2_graph);
    }


    gStyle->SetOptStat(0);
    TCanvas *c = new TCanvas("c", "", 1000, 1000);
    TLegend *l = new TLegend(0.75, 0.75, 0.9, 0.9);
    c->SetLeftMargin(0.15);
    TH1 *dummy = new TH1I("charged_particle_v2", "Charged Particle V2", 2, 0, 20);
    dummy->GetYaxis()->SetRangeUser(-0.05, 0.1);
    dummy->GetYaxis()->SetTitle("v_{2}^{obs}");
    dummy->GetXaxis()->SetTitle("Centrality");
    dummy->Draw("");
    for (int i = 0; i < v2_data->size(); i++) {
        (*v2_data)[i]->SetLineColor(colors[i]);
        (*v2_data)[i]->Draw("l same");
        l->AddEntry((*v2_data)[i], Form("Momentum: %.1f-%.1f", pt_bin_width * i * v2_hist->GetYaxis()->GetBinWidth(1), pt_bin_width * (i+1) * v2_hist->GetYaxis()->GetBinWidth(1)));
    }
    l->Draw();
    c->SaveAs("plots/particle_v2.png");

}