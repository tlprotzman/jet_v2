/*
Tristan Protzman
Lehigh University
February 22, 2022
tlprotzman@gmail.com

Produces (hopefully) pretty QA plots from the jet v2 analysis
*/

#include <TROOT.h>
#include <TTree.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH1D.h>
#include <THStack.h>
#include <TLegend.h>
#include <TFile.h>
#include <TMath.h>
#include <TChain.h>
#include <TStyle.h>
#include <TText.h>
#include <TLatex.h>
#include <TH2.h>
#include <TH2D.h>
#include <TCanvas.h>

#include "setup.h"
#include "histogram_data.h"
#include "histogram_package.h"

#include <iostream>

int main(int argc, char **argv) {
    TFile *file = TFile::Open("out.root");
    TTree *tree = file->Get<TTree>("jet_data");
    jet_tree_data datum;
    read_tree(tree, &datum);

    TLatex *latex = new TLatex;
    latex->SetTextSize(0.03);

    // Histograms
    // Jet Spectra, all jets
    histogram_package jet_pt_spectra;
    jet_pt_spectra.set_save_location("plots/nocut_jet_pt_spectrum");
    jet_pt_spectra.set_log_y();
    jet_pt_spectra.set_x_title("P_{t}");
    jet_pt_spectra.set_y_title("#frac{dN}{dP_{t}}");

    jet_pt_spectra.add_histogram(new histogram_data(file->Get<TH1>("jet_momentum"), "Jet Mometum", kBlue));
    jet_pt_spectra.add_histogram(new histogram_data(file->Get<TH1>("jet_subtracted_momentum"), "Jet Mometum, Background Subtracted", kRed));

    jet_pt_spectra.set_legend(true, 0.48, 0.7);
    jet_pt_spectra.set_x_range(-15, 40);
    jet_pt_spectra.add_drawable(latex->DrawLatexNDC(0.48, 0.83, "Jet Momentum Spectra"));
    
    jet_pt_spectra.draw();
    

    // Track Spectra
    histogram_package track_pt_spectra;
    track_pt_spectra.set_save_location("plots/nocut_track_pt_spectra");
    track_pt_spectra.set_log_y();
    track_pt_spectra.set_x_title("P_{t}");
    track_pt_spectra.set_y_title("#frac{dN}{dP_{t}}");

    track_pt_spectra.add_histogram(new histogram_data(file->Get<TH1>("track_momentum"), "Jet Mometum", kBlue));
    track_pt_spectra.set_x_range(-15, 40);
    track_pt_spectra.add_drawable(latex->DrawLatexNDC(0.48, 0.83, "Track Momentum Spectra"));
    
    track_pt_spectra.draw();

    // Location Spectra
    TH2 *jet_loc = file->Get<TH2>("jet_loc");
    TH1 *jet_eta = jet_loc->ProjectionX();
    TH1 *jet_phi = jet_loc->ProjectionY();

    TH2 *track_loc = file->Get<TH2>("track_loc");
    TH1 *track_eta = track_loc->ProjectionX();
    TH1 *track_phi = track_loc->ProjectionY();

    histogram_package eta_spectra, phi_spectra;
    eta_spectra.set_save_location("plots/nocuts_eta_spectra");
    eta_spectra.set_x_title("#eta");
    eta_spectra.set_y_title("#frac{dN}{d#eta}");

    phi_spectra.set_save_location("plots/nocut_phi_spectra");
    phi_spectra.set_x_title("#phi");
    phi_spectra.set_y_title("#frac{dN}{d#phi}");

    eta_spectra.add_histogram(new histogram_data(track_eta, "Tracks", kRed));
    eta_spectra.add_histogram(new histogram_data(jet_eta, "Jets", kBlue));

    phi_spectra.add_histogram(new histogram_data(track_phi, "Tracks", kRed));
    phi_spectra.add_histogram(new histogram_data(jet_phi, "Jets", kBlue));

    eta_spectra.set_legend();
    phi_spectra.set_legend();

    eta_spectra.add_drawable(latex->DrawLatexNDC(0.5, 0.2, "#eta Distribution"));
    phi_spectra.add_drawable(latex->DrawLatexNDC(0.5, 0.2, "#phi Distribution"));

    eta_spectra.draw();
    phi_spectra.draw();

    TCanvas *c = new TCanvas("", "", 1000, 1000);
    gStyle->SetOptStat(0);
    jet_loc->SetTitle("Jet Location");
    jet_loc->SetXTitle("#eta");
    jet_loc->SetYTitle("#phi");
    jet_loc->Draw("colz");
    c->SaveAs("plots/jet_loc.png");
    c->SaveAs("plots/jet_loc.c");
    delete c;

    c = new TCanvas("", "", 1000, 1000);
    gStyle->SetOptStat(0);
    track_loc->SetTitle("Track Location");
    track_loc->SetXTitle("#eta");
    track_loc->SetYTitle("#phi");
    track_loc->Draw("colz");
    c->SaveAs("plots/track_loc.png");
    c->SaveAs("plots/track_loc.c");
    delete c;

    // Event Plane
    histogram_package event_plane;
    event_plane.set_save_location("plots/event_plane");
    event_plane.set_x_title("#phi");
    event_plane.set_y_title("#frac{dN}{d#phi}");

    histogram_data *east_uncorrected = new histogram_data(file->Get<TH1>("east_uncorrected"), "East, uncorrected", kRed);
    double ep_max = east_uncorrected->get_hist()->GetMaximum() * 1.1;
    double ep_min = east_uncorrected->get_hist()->GetMinimum() * 0.9;
    event_plane.set_y_range(ep_min, ep_max);

    east_uncorrected->get_hist()->SetLineStyle(3);
    east_uncorrected->get_hist()->SetLineWidth(2);
    event_plane.add_histogram(east_uncorrected);
    histogram_data *west_uncorrected = new histogram_data(file->Get<TH1>("west_uncorrected"), "West, uncorrected", kBlue);
    west_uncorrected->get_hist()->SetLineStyle(3);
    west_uncorrected->get_hist()->SetLineWidth(2);
    event_plane.add_histogram(west_uncorrected);
    histogram_data *east_phi_corrected = new histogram_data(file->Get<TH1>("east_phi_corrected"), "East, phi corrected", kRed);
    east_phi_corrected->get_hist()->SetLineStyle(2);
    east_phi_corrected->get_hist()->SetLineWidth(2);
    event_plane.add_histogram(east_phi_corrected);
    histogram_data *west_phi_corrected = new histogram_data(file->Get<TH1>("west_phi_corrected"), "West, phi corrected", kBlue);
    west_phi_corrected->get_hist()->SetLineStyle(2);
    west_phi_corrected->get_hist()->SetLineWidth(2);
    event_plane.add_histogram(west_phi_corrected);
    histogram_data *east_phi_psi_corrected = new histogram_data(file->Get<TH1>("east_phi_psi_corrected"), "East, phi and psi corrected", kRed);
    east_phi_psi_corrected->get_hist()->SetLineStyle(1);
    east_phi_psi_corrected->get_hist()->SetLineWidth(2);
    event_plane.add_histogram(east_phi_psi_corrected);
    histogram_data *west_phi_psi_corrected = new histogram_data(file->Get<TH1>("west_phi_psi_corrected"), "West, phi and psi corrected", kBlue);
    west_phi_psi_corrected->get_hist()->SetLineStyle(1);
    west_phi_psi_corrected->get_hist()->SetLineWidth(2);
    event_plane.add_histogram(west_phi_psi_corrected);

    event_plane.set_legend();

    event_plane.set_title("Event Plane");

    event_plane.draw();

}