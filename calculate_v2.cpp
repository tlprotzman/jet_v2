/*
Tristan Protzman
Lehigh University
Feb 9, 2022
tlprotzman@gmail.com
*/

#include "TROOT.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THStack.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TNamed.h"
#include "TStyle.h"
#include "TTree.h"

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
    double event_plane_east, event_plane_west, event_plane_full;
    double jet_eta, jet_phi;
    double jet_momentum_corrected;
    int centrality_class;

    event_tree->SetBranchAddress("event_plane_east", &event_plane_east);
    event_tree->SetBranchAddress("event_plane_west", &event_plane_west);
    event_tree->SetBranchAddress("event_plane_full", &event_plane_full);
    event_tree->SetBranchAddress("jet_eta", &jet_eta);
    event_tree->SetBranchAddress("jet_phi", &jet_phi);
    event_tree->SetBranchAddress("jet_momentum_median_subtracted", &jet_momentum_corrected);
    event_tree->SetBranchAddress("centrality16", &centrality_class);

    // Set up histograms
    TH1D *phi_spectra = new TH1D("phi_spectra", "Phi Spectra", 100, 0, 3.5);
    TH1D *phi_relative = new TH1D("phi_relative", "Phi relative to ep", 15, 0, TMath::PiOver2());
    TH2D *phi_ep = new TH2D("phi_ep", "Phi vs EP", 100, 0, 3.5, 100, 0, 3.5);
    TH1D *momentum_spectra = new TH1D("momentum_spectra", "Momentum Spectra", 30, 0, 40);
    TH1D *event_plane[6] = {new TH1D("event_plane_east", "event_plane_east", 30, 0, TMath::Pi()),
                            new TH1D("event_plane_west", "event_plane_west", 30, 0, TMath::Pi()),
                            new TH1D("event_plane_east_corrected", "event_plane_east_corrected", 30, 0, TMath::Pi()),
                            new TH1D("event_plane_west_corrected", "event_plane_west_corrected", 30, 0, TMath::Pi()),
                            new TH1D("event_plane_average", "event_plane_average", 30, 0, TMath::Pi()),
                            new TH1D("event_plane_difference", "event_plane_difference", 30, -1 * TMath::Pi(), TMath::Pi())};

    TH2D *ep = new TH2D("ep", "ep", 30, 0, TMath::Pi(), 30, -0, TMath::Pi());
    TH2D *ep_corrected = new TH2D("ep_corrected", "ep_corrected", 30, 0, TMath::Pi(), 30, -0, TMath::Pi());

    double calculated_v2 = 0;
    Long64_t events = 0;
    for (Long64_t n = 0; n < event_tree->GetEntries(); n++) {
        event_tree->GetEvent(n);
        if (centrality_class <= 4 || centrality_class >= 12) {  // accept between 40-20% central
            continue;
        }
        
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

        // std::cout << jet_phi << std::endl;
        if (jet_phi > TMath::Pi()) {
            jet_phi -= TMath::Pi();
        }

        phi_ep->Fill(jet_phi, event_plane_full);

        // jet_phi -= TMath::Pi() / 2.;
        phi_spectra->Fill(jet_phi);

        double relative = abs(jet_phi - event_plane_full);
        while (relative > TMath::PiOver2()) {
            relative -= TMath::PiOver2();
        }
        while (relative < -1 * TMath::PiOver2()) {
            relative += TMath::PiOver2();
        }
        // std::cout << relative << std::endl;
        // if (relative > TMath::PiOver2()) {
        //     relative -= TMath::PiOver2();
        // }
        phi_relative->Fill(relative);
        calculated_v2 += cos(2 * relative);
        events++;

        momentum_spectra->Fill(jet_momentum_corrected);
    }
    std::cout << "Calculated V2: " << calculated_v2 / events << std::endl;

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

    // Fitting histogram
    TF1 *v2_fit = new TF1("v2", "[0] * (1 + [1] * 2 * cos(2*x))", 0, TMath::PiOver2());
    v2_fit->SetParameters(1, 0);
    v2_fit->SetParNames("offset", "v2");
    // phi_relative->Scale(1 / phi_relative->GetBinContent(1));
    phi_relative->Fit("v2");

    // Drawing histograms, probably move later?
    draw_anything(phi_spectra, "plots/phi_spectra.png");
    draw_anything(event_plane[0], "plots/ep_east.png");
    draw_anything(event_plane[1], "plots/ep_west.png");
    event_plane[2]->Scale(1 / event_plane[2]->GetEntries());
    event_plane[3]->Scale(1 / event_plane[3]->GetEntries());
    event_plane[2]->SetMinimum(0);
    event_plane[3]->SetMinimum(0);
    draw_anything(event_plane[2], "plots/ep_east_corrected.png", "hist norm");
    draw_anything(event_plane[3], "plots/ep_west_corrected.png", "hist norm");
    draw_anything(event_plane[4], "plots/ep_average.png");
    draw_anything(event_plane[5], "plots/ep_diff.png");
    draw_anything(momentum_spectra, "plots/momentum_spectra.png");
    gStyle->SetOptStat(0);
    draw_anything(ep, "plots/EP.png", "colz");
    draw_anything(phi_ep, "plots/phi_ep.png", "colz");
    ep_corrected->Scale(1 / ep_corrected->GetEntries());
    draw_anything(ep_corrected, "plots/EP_Corrected.png", "colz");
    // return;

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
    c.Print("plots/stacked.png");


    gPad->SetMargin(0.14, 0.9, 0.9, 0.9);
    phi_relative->Draw("e");
    phi_relative->SetXTitle("#Delta#phi=|#phi^{jet}-#Psi_{2}|");
    phi_relative->SetYTitle("#frac{dN}{d#Delta#phi}");
    phi_relative->GetFunction("v2")->DrawF1(-1 * TMath::PiOver2(), TMath::PiOver2(), "C same");
    v2_fit->SetParameters(phi_relative->GetFunction("v2")->GetParameter("offset"), calculated_v2 / events);
    v2_fit->SetLineColor(kBlue);
    v2_fit->DrawF1(-1 * TMath::PiOver2(), TMath::PiOver2(), "same");
    TLatex formula = TLatex();
    formula.SetTextSize(0.03);
    formula.DrawLatexNDC(0.45, 0.85, Form("Equation: Scale * (1 + 2v_{2} cos(2#phi))", phi_relative->GetFunction("v2")->GetParameter("offset"), phi_relative->GetFunction("v2")->GetParameter("v2")));
    formula.DrawLatexNDC(0.45, 0.8, Form("red: %.1f(1 + 2 * %.3f cos(2#phi))", phi_relative->GetFunction("v2")->GetParameter("offset"), phi_relative->GetFunction("v2")->GetParameter("v2")));
    formula.DrawLatexNDC(0.45, 0.75, Form("blue: %.1f(1 + 2 * %.3f cos(2#phi))", phi_relative->GetFunction("v2")->GetParameter("offset"), calculated_v2 / events));
    
    formula.DrawLatexNDC(0.2, 0.2, "Red: Fitted to data");
    formula.DrawLatexNDC(0.2, 0.15, "Blue: v_{2}= #LTcos(2#Delta#phi)#GT");
    c.Print("plots/phi_relative.png");
}