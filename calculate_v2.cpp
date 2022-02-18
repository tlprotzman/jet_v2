/*
Tristan Protzman
Lehigh University
Feb 9, 2022
tlprotzman@gmail.com
*/

#include "TROOT.h"
#include "TFile.h"
#include "TNamed.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "TMath.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "THStack.h"

#include <iostream>
#include <string>


void draw_anything(TNamed *hist, std::string name, std::string opt="") {
    TCanvas c("", "", 1000, 1000);
    hist->Draw(opt.c_str());
    c.Print(name.c_str());
}

void calculate_v2(std::string infile_path = "out.root") {
    // open file and get tree
    ROOT::EnableImplicitMT();
    TFile *infile = TFile::Open(infile_path.c_str());
    TTree *event_tree = infile->Get<TTree>("jet_data");


    // Set up needed branches
    double event_plane_east, event_plane_west;
    double jet_eta, jet_phi;
    double jet_momentum_corrected;

    event_tree->SetBranchAddress("event_plane_east", &event_plane_east);
    event_tree->SetBranchAddress("event_plane_west", &event_plane_west);
    // event_tree->SetBranchAddress("jet_eta", &jet_eta);
    // event_tree->SetBranchAddress("jet_phi", &jet_phi);
    // event_tree->SetBranchAddress("jet_momentum_median_subtracted", &jet_momentum_corrected);

    // Set up histograms
    TH1D *phi_spectra = new TH1D("phi_spectra", "Phi Spectra", 100, -1, 7);
    TH1D *phi_relative = new TH1D("phi_relative", "Phi relative to ep", 30, 0, TMath::TwoPi());
    TH1D *event_plane[6] = {new TH1D("event_plane_east", "event_plane_east", 30, 0, TMath::TwoPi()),
                            new TH1D("event_plane_west", "event_plane_west", 30, 0, TMath::TwoPi()),
                            new TH1D("event_plane_east_corrected", "event_plane_east_corrected", 30, 0, TMath::TwoPi()),
                            new TH1D("event_plane_west_corrected", "event_plane_west_corrected", 30, 0, TMath::TwoPi()),
                            new TH1D("event_plane_average", "event_plane_average", 30, 0, TMath::TwoPi()),
                            new TH1D("event_plane_difference", "event_plane_difference", 30, -1 * TMath::Pi(), TMath::Pi())};

    TH2D *ep = new TH2D("ep", "ep", 30, 0, TMath::TwoPi(), 30, -0, TMath::TwoPi());
    TH2D *ep_corrected = new TH2D("ep_corrected", "ep_corrected", 30, 0, TMath::TwoPi(), 30, -0, TMath::TwoPi());

    for (Long64_t n = 0; n < event_tree->GetEntries(); n++) {
        event_tree->GetEvent(n);
        
        float average_ep = 0.5 * (event_plane_east + event_plane_west);

        // Fill ep histograms
        event_plane[0]->Fill(event_plane_east);
        event_plane[1]->Fill(event_plane_west);
        event_plane[4]->Fill(average_ep);
        double diff = event_plane_east - event_plane_west;
        if (diff < 0) {
            diff += TMath::TwoPi();
        } 
        event_plane[5]->Fill(event_plane_east - event_plane_west);
        ep->Fill(event_plane_east, event_plane_west);

        phi_spectra->Fill(jet_phi);

        double relative = jet_phi - event_plane_west;
        if (relative < 0) {
            relative += TMath::TwoPi();
        }
        phi_relative->Fill(relative);
    }

    // Apply corrections to EPD to account for bias from missing TPC sector
    double **correction_factor = (double**) malloc(2 * sizeof(double*));
    correction_factor[0] = (double*) malloc((event_plane[0]->GetNbinsX() + 1 )* sizeof(double));
    correction_factor[1] = (double*) malloc((event_plane[1]->GetNbinsX() + 1 )* sizeof(double));  // will break if not equal number of bins, whoops
    
    for (size_t i = 0; i < event_plane[0]->GetNbinsX() + 1; i++) {
        correction_factor[0][i] = event_plane[0]->GetBinContent(i) / event_plane[0]->GetEntries();
        correction_factor[1][i] = event_plane[1]->GetBinContent(i) / event_plane[1]->GetEntries();
        std::cout << correction_factor[0][i] << "\t" << correction_factor[1][i] << std::endl;
    }

    for (Long64_t n = 0; n < event_tree->GetEntries(); n++) {
        event_tree->GetEvent(n);
        int east_bin = 1;
        int west_bin = 1;
        while (event_plane_east > event_plane[0]->GetBinLowEdge(east_bin + 1)) {east_bin++;}
        while (event_plane_west > event_plane[1]->GetBinLowEdge(west_bin + 1)) {west_bin++;}
        event_plane[2]->Fill(event_plane_east, 1 / correction_factor[0][east_bin]);
        // std::cout << "weight: " << 1 / correction_factor[0][east_bin] << "\tI: " << east_bin << std::endl;
        event_plane[3]->Fill(event_plane_west, 1 / correction_factor[1][west_bin]);
        ep_corrected->Fill(event_plane_east, event_plane_west, (1 / correction_factor[0][east_bin]) * (1 / correction_factor[1][west_bin]));
    }


    // Drawing histograms, probably move later?
    // draw_anything(phi_spectra, "phi_spectra.png");
    // draw_anything(phi_relative, "phi_relative.png");
    draw_anything(event_plane[0], "ep_east.png");
    draw_anything(event_plane[1], "ep_west.png");
    event_plane[2]->Scale(1 / event_plane[2]->GetEntries());
    event_plane[3]->Scale(1 / event_plane[3]->GetEntries());
    event_plane[2]->SetMinimum(0);
    event_plane[3]->SetMinimum(0);
    draw_anything(event_plane[2], "ep_east_corrected.png", "hist norm");
    draw_anything(event_plane[3], "ep_west_corrected.png", "hist norm");
    draw_anything(event_plane[4], "ep_average.png");
    draw_anything(event_plane[5], "ep_diff.png");
    gStyle->SetOptStat(0);
    draw_anything(ep, "EP.png", "colz");
    ep_corrected->Scale(1 / ep_corrected->GetEntries());
    draw_anything(ep_corrected, "EP_Corrected.png", "colz");
    return;

    TCanvas c("", "", 1000, 1000);
    
    THStack *stack = new THStack();
    TLegend *legend = new TLegend();
    phi_spectra->SetLineColor(kBlue);
    stack->Add(phi_spectra);
    legend->AddEntry(phi_spectra, "Spectra");
    phi_relative->SetLineColor(kRed);
    stack->Add(phi_relative);
    legend->AddEntry(phi_relative, "Relative to Event Plane");


    stack->Draw("nostack E");
    legend->Draw();
    c.Print("stacked.png");

}