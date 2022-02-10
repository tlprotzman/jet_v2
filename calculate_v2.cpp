/*
Tristan Protzman
Lehigh University
Feb 9, 2022
tlprotzman@gmail.com
*/

#include "TROOT.h"
#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"
#include "TMath.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "THStack.h"

#include <iostream>
#include <string>


void draw_hist(TH1 *hist, std::string name) {
    TCanvas c("", "", 1000, 1000);
    hist->Draw();
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
    event_tree->SetBranchAddress("jet_eta", &jet_eta);
    event_tree->SetBranchAddress("jet_phi", &jet_phi);
    event_tree->SetBranchAddress("jet_momentum_corrected", &jet_momentum_corrected);

    // Set up histograms
    TH1D *phi_spectra = new TH1D("phi_spectra", "Phi Spectra", 30, 0, TMath::TwoPi());
    TH1D *phi_relative = new TH1D("phi_relative", "Phi relative to ep", 30, 0, TMath::TwoPi());
    TH1D *event_plane[4] = {new TH1D("event_plane_east", "event_plane_east", 30, 0, TMath::TwoPi()),
                            new TH1D("event_plane_west", "event_plane_west", 30, 0, TMath::TwoPi()),
                            new TH1D("event_plane_average", "event_plane_average", 30, 0, TMath::TwoPi()),
                            new TH1D("event_plane_difference", "event_plane_difference", 30, -1 * TMath::Pi(), TMath::Pi())};

    for (Long64_t n = 0; n < event_tree->GetEntries(); n++) {
        event_tree->GetEvent(n);
        
        float average_ep = 0.5 * (event_plane_east + event_plane_west);

        // Fill ep histograms
        event_plane[0]->Fill(event_plane_east);
        event_plane[1]->Fill(event_plane_west);
        event_plane[2]->Fill(average_ep);
        event_plane[3]->Fill(event_plane_east - event_plane_west);

        phi_spectra->Fill(jet_phi);

        double relative = jet_phi - event_plane_east;
        if (relative < 0) {
            relative += TMath::TwoPi();
        }
        phi_relative->Fill(relative);
    }

    // Drawing histograms, probably move later?
    draw_hist(phi_spectra, "phi_spectra.png");
    draw_hist(phi_relative, "phi_relative.png");
    draw_hist(event_plane[0], "ep_east.png");
    draw_hist(event_plane[1], "ep_west.png");
    draw_hist(event_plane[2], "ep_average.png");
    draw_hist(event_plane[3], "ep_diff.png");

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