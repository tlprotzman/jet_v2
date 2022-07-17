#include "particle_anisotropy.h"

#include <TROOT.h>
#include <TH1.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLatex.h>

#include "fastjet/PseudoJet.hh"

#include <vector>

int JET_NUMBER = 0;

double JET_R = 0.2; // This shouldn't be hardcoded...

particle_anisotropy::particle_anisotropy() {
    return;
}

particle_anisotropy::~particle_anisotropy() {
    return;
}

void particle_anisotropy::calculate_v2(std::vector<fastjet::PseudoJet> &tracks, fastjet::PseudoJet &leading_jet, double event_plane, particle_anisotropy &flow_parameters, int centrality) {
    double j_eta = leading_jet.eta();
    double j_phi = leading_jet.phi();

    // Create list of tracks outside leading jet eta range
    std::vector<fastjet::PseudoJet> *candidate_tracks = new std::vector<fastjet::PseudoJet>();
    for (fastjet::PseudoJet track : tracks) {
        double d_eta = abs(track.eta() - j_eta);
        if (d_eta > JET_R) {
            candidate_tracks->push_back(track);
        }
    }
    std::cout << "Calculating flow quantities with " << candidate_tracks->size() << " tracks\n";

    // Populate azimuthal angle distribution
    TH1 *particle_distribution = new TH1D("particle_distribution", "", sqrt(candidate_tracks->size()), 0, TMath::TwoPi());
    for (fastjet::PseudoJet track : tracks) {
        double relative = track.phi() - event_plane;
        if (relative < 0) {
            relative += TMath::TwoPi();
        }
        particle_distribution->Fill(relative, track.pt());  // Sum of momentum in azimuthal bins
    }

    // Extract flow parameters
    TF1 *flow_description = new TF1("flow_description", "[0] * (1 + 2 * ([1] * cos(2 * x)))", 0, TMath::TwoPi(), "Q");
    flow_description->SetParNames("rho_0", "v_2");
    
    particle_distribution->Fit(flow_description, "q");
    flow_parameters.rho = flow_description->GetParameter("rho_0");
    flow_parameters.v2 = flow_description->GetParameter("v_2");
    flow_parameters.v3 = 0;//flow_description->GetParameter("v_3");

    if (centrality == 0) {
        this->print_particle_distribution(particle_distribution);
        JET_NUMBER++;
    }

    delete flow_description;
    delete particle_distribution;
    delete candidate_tracks;

    return;
}

void particle_anisotropy::print_particle_distribution(TH1 *hist) {
    TCanvas *c = new TCanvas("", "", 1000, 1000);
    hist->Draw("e");
    TF1 *fit = hist->GetFunction("flow_description");
    fit->Draw("c same");
    TF1 *constant = new TF1("constant_term", Form("%f", fit->GetParameter("rho_0")), 0, TMath::TwoPi());
    constant->SetLineColor(kMagenta);
    constant->Draw("c same");
    TF1 *v2 = new TF1("v2_term", Form("%f * (1 + 2 * %f * cos(2 * x))", fit->GetParameter("rho_0"), fit->GetParameter("v_2")), 0, TMath::TwoPi());
    v2->SetLineColor(kBlue);
    v2->Draw("c same");
    // TF1 *v3 = new TF1("v3_term", Form("%f * (1 + 2 * %f * cos(3 * x))", fit->GetParameter("rho_0"), fit->GetParameter("v_3")), 0, TMath::TwoPi());
    // v3->SetLineColor(kGreen+3);
    // v3->Draw("c same");

    fit->GetXaxis()->SetTitle("#Delta#phi");
    fit->GetYaxis()->SetTitle("#rho_{charged} (GeV/c)");

    TLatex *description = new TLatex(0.15, 0.2, Form("#rho_{0}: %.5f     v_{2}: %.5f     v_{3}: %.5f", fit->GetParameter("rho_0"), fit->GetParameter("v_2"), 0/*fit->GetParameter("v_3")*/));
    description->SetTextSize(0.03);
    description->SetNDC();

    description->Draw();
    c->SaveAs(Form("plots/event_%i.png", JET_NUMBER));

    delete constant;
    delete v2;
    // delete v3;
    delete description;

    delete c;
}