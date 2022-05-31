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
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"

#include "jet_tree.h"
#include "event_tree.h"

#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <vector>

// const int num_cent_bins = 16;
// char *cent_bins[num_cent_bins] = {"75-80%", "70-75%", "65-70%", "60-65%", "55-60%", "50-55%", "45-50%", "40-45%", "35-40%", "30-35%", "25-30%", "20-25%", "15-20%", "10-15%", "5-10%", "0-5%"};
int colors[16] = {kAzure + 10, kAzure,
                         kViolet + 10, kViolet-3, kViolet - 9,
                         kPink + 10, kPink,
                         kOrange + 10, kOrange, kOrange - 9,
                         kSpring + 10, kSpring,
                         kTeal + 10, kTeal, kTeal - 9, kCyan + 4
                         };
// const int num_cent_bins = 4;
// char *cent_bins[num_cent_bins] = {"60-80%", "40-60%", "20-40%", "0-20%"};
const int num_cent_bins = 8;
char *cent_bins[num_cent_bins] = {"70-80%", "60-70%", "50-60%", "40-50%", "30-40%", "20-30%", "10-20%", "0-10%",};
// int colors[num_cent_bins] = {kRed, kOrange, kBlue, kViolet};

const int num_pt_bins = 4;

// Forward declarations
TChain* load_files(std::string file_list);
TH3* setup_histogram();
void event_loop(TChain *chain, Event_Tree *event_tree, Jet_Tree *jet_tree, Jet_Tree *hardcore_jet_tree, TH3 *relative_angle, double *ep_resolution);
void save_histogram(TH3 *hist, std::string outfile_name);
int calculate_v2(TH1 *dphi, int cent_bin, int pt_bin);
double*** bin_loop(TH3 *relative, double *ep_resolution);
void draw_plot(TH1 *dphi, int cent_bin, int pt_bin);
void plot_v2_values(double ***v2);
void plot_ep_reso(double *ep_reso);

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
    double *ep_resolution = (double*)malloc(num_cent_bins * sizeof(double));
    for (int i = 0; i < num_cent_bins; i++) {
        ep_resolution[i] = 0;
    }
    event_loop(chain, event_tree, jet_tree, hardcore_jet_tree, relative_angle, ep_resolution);

    double ***v2 = bin_loop(relative_angle, ep_resolution);
    plot_v2_values(v2);
    plot_ep_reso(ep_resolution);


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

    int n_pt_bins = num_pt_bins;
    double pt_low = 0;
    double pt_high = 10 * n_pt_bins; // 5 GeV bins

    int n_centrality_bins = num_cent_bins;
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

bool hardcore_matched(Jet_Tree *jet_tree, Jet_Tree *hardcore_jet_tree, int index) {     // Currently hardcore jets can be reused...
    double j_phi = jet_tree->jet_phi[index];
    double j_eta = jet_tree->jet_eta[index];
    for (uint i = 0; i < hardcore_jet_tree->num_jets; i++) {
        double d_phi = abs(j_phi - hardcore_jet_tree->jet_phi[i]);
        double d_eta = abs(j_eta - hardcore_jet_tree->jet_eta[i]);
        double dr = sqrt(d_phi * d_phi + d_eta * d_eta);
        if (dr < 0.3) {
            // std::cout << "matched!\n";
            return true;
        }
    }
    // std::cout << "not matched :(\n";
    return false;
}

void event_loop(TChain *chain, Event_Tree *event_tree, Jet_Tree *jet_tree, Jet_Tree *hardcore_jet_tree, TH3 *relative_angle, double *ep_resolution) {
    int *n_events = (int*) malloc(num_cent_bins * sizeof(int));
    for (int i = 0; i < num_cent_bins; i++) {
        n_events[i] = 0;
    }

    for (Long64_t n = 0; n < chain->GetEntries(); n++) {
        // if (n > 1000000) {
        //     return;
        // }
        if (n % 100000 == 0) {
            printf("Processed %.1fM events\n", n / 1000000.);
        }
        chain->GetEvent(n); // Get the next events

        if (event_tree->centrality >= 16 || event_tree->centrality < 0) {
            continue; // Either very peripheral or can't determine 
        }
        n_events[event_tree->centrality]++;
        ep_resolution[event_tree->centrality] += cos(2 * (event_tree->ep_east - event_tree->ep_west));

        // For each jet in the event...
        for (int i = 0; i < jet_tree->num_jets; i++) {
            if (!hardcore_matched(jet_tree, hardcore_jet_tree, i)) {
                continue;
            }
            if (jet_tree->jet_charged_z[i] > 0.95) {    // Skip jets mostly composed of single track
                continue;
            }
            if (abs(jet_tree->jet_eta[i]) > 0.7) {      // Skip jets within R of detector edge
                continue;
            }

            double relative = jet_tree->jet_phi[i] - event_tree->ep_east;
            if (relative < 0) {
                relative += TMath::TwoPi();
            }
            if (relative > TMath::Pi()) {
                relative -= TMath::Pi();
            }
            double background = jet_tree->rho;
            double assumed_background_v2 = 0.055;
            background = background * (1 + assumed_background_v2 * cos(2 * (relative)));
            relative_angle->Fill(relative, jet_tree->jet_pt[i] - background * jet_tree->jet_area_pt[i], event_tree->centrality);
        }
    }

    for (int i = 0; i < num_cent_bins; i++) {
        ep_resolution[i] /= n_events[i];
        ep_resolution[i] = sqrt(ep_resolution[i]) * sqrt(2);
        std::cout << "ep reso for " << cent_bins[i] << ": " << ep_resolution[i] << " with " << n_events[i] << "events" << std::endl;
    }
}

void save_histogram(TH3 *hist, std::string outfile_name) {
    TFile *outfile = new TFile(outfile_name.c_str(), "RECREATE");
    hist->SetDirectory(outfile);
    hist->Write();
    outfile->Close();
    delete outfile;
}

// Fits the data and extracts v2, returns 0 if no entries in bin
int calculate_v2(TH1 *dphi, int cent_bin, int pt_bin) {
    if (dphi->GetEntries() == 0) {
        return 0;
    }
    TF1 *v2_fit = new TF1(Form("v2_bin_%d_%d", cent_bin, pt_bin), "[0] * (1 + [1] * 2 * cos(2*x))", 0, TMath::Pi()); // Fit distribution on [0, pi]
    v2_fit->SetParameters(1, 0);
    v2_fit->SetParNames("offset", "v2");
    dphi->Fit(Form("v2_bin_%d_%d", cent_bin, pt_bin)); // Run fitting
    draw_plot(dphi, cent_bin, pt_bin);
    return 1;
}

// Loops over and determins v2 for each bin
double*** bin_loop(TH3 *relative, double *ep_resolution) {
    //v2[0][:][:] are v2 values, v2[1][:][:] are the associated errors
    double ***v2 = (double***)malloc(2 * sizeof(double**)); 
    v2[0] = (double**)malloc(num_cent_bins * sizeof(double*)); //v2[:][cent_bin][:]
    v2[1] = (double**)malloc(num_cent_bins * sizeof(double*)); //v2[:][cent_bin][:]
    for (int cent_bin = 1; cent_bin <= num_cent_bins; cent_bin++) {
        v2[0][cent_bin - 1] = (double*)malloc(8 * sizeof(double)); //v2[:][:][pt_bin]
        v2[1][cent_bin - 1] = (double*)malloc(8 * sizeof(double)); //v2[:][:][pt_bin]

        relative->GetZaxis()->SetRange(cent_bin, cent_bin); // Restrict centrality range
        for (int pt_bin = 1; pt_bin <= relative->GetYaxis()->GetNbins(); pt_bin++) {
            relative->GetYaxis()->SetRange(pt_bin, pt_bin); // Restrict pt range
            TH1 *dphi = relative->Project3D("x");
            if (calculate_v2(dphi, cent_bin, pt_bin)) {
                TF1 *func = dphi->GetFunction(Form("v2_bin_%d_%d", cent_bin, pt_bin));
                v2[0][cent_bin - 1][pt_bin - 1] = func->GetParameter("v2") / ep_resolution[cent_bin - 1];
                v2[1][cent_bin - 1][pt_bin - 1] = func->GetParError(func->GetParNumber("v2")) / ep_resolution[cent_bin - 1];
                // std::cout << "v2: " << << "\tError: " << << std::endl;
            }
            else {
                v2[0][cent_bin - 1][pt_bin - 1] = -999;
                v2[1][cent_bin - 1][pt_bin - 1] = 0;
            }
        }
    }
    return v2;
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

void plot_v2_values(double ***v2) {
    TCanvas *c = new TCanvas("", "", 2000, 1000);
    // c->Divide(2, 1);
    // c->cd(1);
    TLegend *legend1 = new TLegend(0.1, 0.1, 0.25, 0.3);
    TLegend *legend2 = new TLegend(0.1, 0.1, 0.25, 0.3);

    double *x = (double*)malloc(num_pt_bins * sizeof(double));
    double *ex = (double*)malloc(num_pt_bins * sizeof(double));
    for (int i = 0; i < num_pt_bins; i++) {
        x[i] = 10 * i + 4.7;
        ex[i] = 1;
    }

    
    std::vector<TGraphErrors*> v2_graphs;
    bool first = true;
    for (int cent_bin = 0; cent_bin < num_cent_bins; cent_bin ++) {
        std::cout << "look at me!" << v2[1][cent_bin][1] << std::endl;
        TGraphErrors *v2_plot = new TGraphErrors(num_pt_bins, x, v2[0][cent_bin], ex, v2[1][cent_bin]);
        // TGraphErrors *v2_plot = new TGraphErrors(num_pt_bins, x, x, ex, ex);
        v2_plot->SetLineColor(colors[cent_bin]);
        v2_plot->SetMarkerStyle(kDot);
        v2_plot->SetMarkerSize(15);
        v2_plot->SetMarkerColor(colors[cent_bin]);
        v2_plot->SetLineWidth(2);
        // if (cent_bin < 8) {
            legend1->AddEntry(v2_plot, cent_bins[cent_bin]);
        // }
        // else {
        //     if (first) {
        //         first = false;
        //         for (int i = 0; i < 8; i++) {
        //             x[i] = 5 * i + 1.7;
        //         }
        //     }
        //     legend2->AddEntry(v2_plot, cent_bins[cent_bin]);
        // }
        v2_graphs.push_back(v2_plot);
        for (int i = 0; i < num_pt_bins; i++) {
            x[i] += 0.2;
        }
    }

    v2_graphs[0]->Draw("alp");
    v2_graphs[0]->SetTitle("Jet v_{2}, hardcore matched");
    v2_graphs[0]->GetXaxis()->SetTitle("p_{T}");
    v2_graphs[0]->GetYaxis()->SetTitle("v_{2}");
    v2_graphs[0]->GetYaxis()->SetRangeUser(-0.2, 0.4);
    for (int i = 1; i < num_cent_bins; i++) {
        v2_graphs[i]->Draw("pl same");
    }
    legend1->Draw();

    // c->cd(2);
    // v2_graphs[8]->Draw("ap");
    // v2_graphs[8]->SetTitle("Jet v_{2}");
    // v2_graphs[8]->GetXaxis()->SetTitle("p_{T}");
    // v2_graphs[8]->GetYaxis()->SetTitle("v_{2}");
    // v2_graphs[8]->GetYaxis()->SetRangeUser(-0.2, 0.4);
    // for (int i = 9; i < v2_graphs.size(); i++) {
    //     v2_graphs[i]->Draw("p same");

    // }
    // legend2->Draw();
    c->SaveAs("plots/v2_main_hc.png");
    c->SaveAs("plots/v2_main_hc.c");
}

void plot_ep_reso(double *ep_reso) {
    double x[num_cent_bins];
    for (int i = 0; i < num_cent_bins; i++) {
        x[i] = i;
        std::cout << "Bin: " << i << "\tReso: " << ep_reso[i] << std::endl;
    }
    TGraph *ep_reso_graph = new TGraph(num_cent_bins, x, ep_reso);
    ep_reso_graph->SetMarkerSize(1);
    ep_reso_graph->SetMarkerStyle(kFullCircle);
    ep_reso_graph->GetXaxis()->SetNdivisions(100 + (num_cent_bins / 2));
    for (uint32_t i = 0; i < num_cent_bins / 2; i++) {
        ep_reso_graph->GetXaxis()->ChangeLabel(i + 1, 315, -1, -1, -1, -1, cent_bins[2 * i]);
    }
    TCanvas *c = new TCanvas("", "", 1000, 1000);
    ep_reso_graph->Draw("apw");
    gPad->SetGrid(1, 1);
    ep_reso_graph->GetXaxis()->SetLabelOffset(0.04);
    ep_reso_graph->GetXaxis()->SetLabelSize(0.025);

    ep_reso_graph->SetTitle("Event Plane Resolution");
    ep_reso_graph->GetXaxis()->SetTitle("Centrality");
    ep_reso_graph->GetYaxis()->SetTitle("Resolution");

    c->SaveAs("plots/ep_reso.png");
    c->SaveAs("plots/ep_reso.c");

}

void merge_v2_bins(TH3 *relative) {

}