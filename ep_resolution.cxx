/*
Calculates the EP resolution as a function of centrality class
for rapid iteration on weighting methods
Tristan Protzman
Lehigh University
tlp220@lehigh.edu
April 22, 2022
*/

#include "TROOT.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TPad.h"
#include "TStyle.h"

#include <stdint.h>
#include <iostream>

void ep_resolution() {
    // Let's load the histograms we need
    TFile *infile = TFile::Open("out.root");
    TH2D *ep_correlation = infile->Get<TH2D>("ep_correlation");
    TH2D *ep_resolution  = infile->Get<TH2D>("ep_resolution");

    // This is where the calculations will go - I'll get to those later
    int num_centrallity_classes = ep_resolution->GetNbinsY();
    std::cout << num_centrallity_classes << " centrallity classes" << std::endl;
    double y[17];
    for (uint32_t cent_bin = 0; cent_bin < num_centrallity_classes - 1; cent_bin++) { // loop over centrallity bins
        uint32_t n_events = 0;
        double ep_res = 0;
        for (uint32_t diff_bin = 0; diff_bin < ep_resolution->GetNbinsX(); diff_bin++) {
            uint32_t bin_entries = ep_resolution->GetBinContent(diff_bin, cent_bin);
            n_events += bin_entries;
            ep_res += bin_entries * cos(2 * ep_resolution->GetXaxis()->GetBinCenter(diff_bin));
        }
        if (n_events) {
            ep_res /= n_events;
        }
        else {
            ep_res = 0;
        }
        ep_res = sqrt(abs(ep_res)) * sqrt(2);
        std::cout << cent_bin << "\t" << ep_res << std::endl;
        y[cent_bin] = ep_res;
    }


    // Drawing what I come up with
    TCanvas *canvas = new TCanvas("", "", 3000, 1000);
    canvas->Divide(3, 1);
    gStyle->SetOptStat(11); // Maybe?  I can never remember
    
    // Correlation plot
    canvas->cd(1);
    canvas->SetLogz();
    ep_correlation->Draw("colz");

    // Resolution plot - 2D
    canvas->cd(2);
    canvas->SetLogz();
    ep_resolution->Draw("colz");


    // 1D Resolution plot
    canvas->cd(3);
    double x[17] = {15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 16};
    TGraph *ep_resolution_binned = new TGraph(17, x, y);
    ep_resolution_binned->SetMarkerSize(1);
    ep_resolution_binned->SetMarkerStyle(kFullCircle);
    const char * label[17] = {"0-5", "5-10", "10-15", "15-20", "20-25", "25-30", "30-35", "35-40", "40-45", "45-50", "50-55", "55-60", "60-65", "65-70", "70-75", "75-80", "80-100"};
    ep_resolution_binned->Draw("APW");
    // ep_resolution_binned->GetXaxis()->SetNdivisions(17);
    for (uint32_t i = 1; i < 10; i++) {
        ep_resolution_binned->GetXaxis()->ChangeLabel(i, 315, -1, -1, -1, -1, label[(i - 1) * 2]);
    }
    gPad->SetGrid(1, 1);
    ep_resolution_binned->GetXaxis()->SetLabelOffset(0.04);
    ep_resolution_binned->GetXaxis()->SetLabelSize(0.025);

    
    // Saving plot
    canvas->SaveAs("plots/ep_resolution_debugging.png");
    canvas->SaveAs("plots/ep_resolution_debugging.C");

    // Cleanup - lets keep some decent habits
    infile->Close();
    delete canvas;
        
}