#include "TROOT.h"
#include "TChain.h"
#include "TF1.h"
#include "TH1.h"
#include "TH3.h"
#include "TH3D.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"

#include "jet_tree.h"
#include "event_tree.h"

#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>


// Forward declarations
TChain* load_files(std::string file_list);
TH3* setup_histogram();
void event_loop(TChain *chain, Event_Tree *event_tree, Jet_Tree *jet_tree, TH3 *relative_angle);
void save_histogram(TH3 *hist, std::string outfile_name);
void calculate_v2(TH1 *dphi, int cent_bin, int pt_bin);
void bin_loop(TH3 *relative);
void draw_plot(TH1 *dphi, int cent_bin, int pt_bin);

int main(int argc, char **argv) {
    TChain *chain = load_files("in.list");

    

    // Set up trees for reading
    Event_Tree *event_tree = new Event_Tree(chain, "event");    // Event wise quantities
    event_tree->readable_tree();
    Jet_Tree *hardcore_jet_tree = new Jet_Tree(chain, "hc");    // Hardcore jets
    hardcore_jet_tree->readable_tree();
    Jet_Tree *jet_tree = new Jet_Tree(chain, "all");            // All jets
    jet_tree->readable_tree();
    // return 0;

    // Create momentum and centrality bins
    TH3 *relative_angle = setup_histogram();

    // Process events
    event_loop(chain, event_tree, jet_tree, relative_angle);

    bin_loop(relative_angle);

    save_histogram(relative_angle, "relative_angles.root");
    // delete relative_angle;
    // delete event_tree;
    // delete jet_tree;
    // delete hardcore_jet_tree;
    // delete chain;


}

TChain* load_files(std::string file_list) {
    // Load trees to analyze
    TChain *chain = new TChain("jet_data");
    std::ifstream files("in.list");
    std::string filePath;
    while (files >> filePath) {
        chain->AddFile(filePath.c_str());
    }
    files.close();
    return chain;
}

TH3* setup_histogram() {
    int n_phi_bins = 50;
    double phi_low = 0;
    double phi_high = TMath::Pi();

    int n_pt_bins = 8;
    double pt_low = 0;
    double pt_high = 5 * n_pt_bins; // 5 GeV bins

    int n_centrality_bins = 16;
    double centrality_low = -0.5;
    double centrality_high = 15.5;

    TH3 *h = new TH3D("relative", "relative",
                      n_phi_bins, phi_low, phi_high,
                      n_pt_bins, pt_low, pt_high,
                      n_centrality_bins, centrality_low, centrality_high);

    h->GetXaxis()->SetTitle("Phi");
    h->GetYaxis()->SetTitle("Momentum");
    h->GetZaxis()->SetTitle("Centrality");

    return h;
}

void event_loop(TChain *chain, Event_Tree *event_tree, Jet_Tree *jet_tree, TH3 *relative_angle) {
    for (Long64_t n = 0; n < chain->GetEntries(); n++) {
        if (n % 100000 == 0) {
            printf("Processed %.1fM events\n", n / 1000000.);
        }
        chain->GetEvent(n); // Get the next events

        if (event_tree->centrality >= 16 || event_tree->centrality < 0) {
            continue; // Either very peripheral or can't determine 
        }
        
        // For each jet in the event...
        for (int i = 0; i < jet_tree->num_jets; i++) {
            double relative = jet_tree->jet_phi[i] - event_tree->ep_east;
            if (relative < 0) {
                relative += TMath::TwoPi();
            }
            if (relative > TMath::Pi()) {
                relative -= TMath::Pi();
            }
            relative_angle->Fill(relative, jet_tree->jet_pt_median_subtracted[i], event_tree->centrality);
        }
    }
}

void save_histogram(TH3 *hist, std::string outfile_name) {
    TFile *outfile = new TFile(outfile_name.c_str(), "RECREATE");
    hist->SetDirectory(outfile);
    hist->Write();
    outfile->Close();
    delete outfile;
}

// Fits the data and extracts v2
void calculate_v2(TH1 *dphi, int cent_bin, int pt_bin) {
    if (dphi->GetEntries() == 0) {
        return;
    }
    TF1 *v2_fit = new TF1(Form("v2_bin_%d_%d", cent_bin, pt_bin), "[0] * (1 + [1] * 2 * cos(2*x))", 0, TMath::Pi()); // Fit distribution on [0, pi]
    v2_fit->SetParameters(1, 0);
    v2_fit->SetParNames("offset", "v2");
    dphi->Fit(Form("v2_bin_%d_%d", cent_bin, pt_bin)); // Run fitting
    draw_plot(dphi, cent_bin, pt_bin);
}

// Loops over and determins v2 for each bin
void bin_loop(TH3 *relative) {
    for (int cent_bin = 1; cent_bin <= relative->GetZaxis()->GetNbins(); cent_bin++) {
        relative->GetZaxis()->SetRange(cent_bin, cent_bin);
        for (int pt_bin = 1; pt_bin <= relative->GetYaxis()->GetNbins(); pt_bin++) {
            relative->GetYaxis()->SetRange(pt_bin, pt_bin);
            TH1 *dphi = relative->Project3D("x");
            calculate_v2(dphi, cent_bin, pt_bin);
        }
    }
}

void draw_plot(TH1 *dphi, int cent_bin, int pt_bin) {
    TCanvas c("", "", 1000, 1000);
    dphi->SetLineColor(kRed);
    gPad->SetMargin(0.14, 0.9, 0.9, 0.9);
    dphi->Draw("e");
    dphi->SetXTitle("#Delta#phi=|#phi^{jet}-#Psi_{2}|");
    dphi->SetYTitle("#frac{dN}{d#Delta#phi}");
    TF1 *func = dphi->GetFunction(Form("v2_bin_%d_%d", cent_bin, pt_bin));
    func->DrawF1(0, TMath::Pi(), "C same");
    TLatex formula = TLatex();
    formula.SetTextSize(0.02);
    // formula.DrawLatexNDC(0.35, 0.85, Form("Equation: Scale * (1 + 2v_{2} cos(2#phi))", dphi->GetFunction("v2")->GetParameter("offset"), dphi->GetFunction("v2")->GetParameter("v2")));
    formula.DrawLatexNDC(0.35, 0.85, Form("%.1f(1 + 2 * %.5f cos(2#Delta#phi))", func->GetParameter("offset"), func->GetParameter("v2")));
    // formula.DrawLatexNDC(0.35, 0.79, "30\%-50\% central");
    // formula.DrawLatexNDC(0.35, 0.79, "30\%-50\% central");
    
    c.Print(Form("plots/slices/png/jet_relative_phi_cent_%d_pt_%d.png", cent_bin, pt_bin));
    c.Print(Form("plots/slices/c/jet_relative_phi_cent_%d_pt_%d.c", cent_bin, pt_bin));

}
