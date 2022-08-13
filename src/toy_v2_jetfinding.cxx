#include <TROOT.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <THStack.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TMath.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TSystem.h>

#include <fastjet/AreaDefinition.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>

#include <iostream>
#include <vector>

double boltzmann_distribution(TRandom *rng, double temperature, double max) {
    double track_pt;
    double p;
    do {
        track_pt = rng->Rndm() * max;    // 10 GeV max track 
    } while (rng->Rndm() > exp(-1 * track_pt / temperature));
    return track_pt;
}


void simulation(int job_id) {

    // Set up model parameters
    int n = 250000; // Number of events to throw
    
    double max_eta = 1;  // Bounds on eta
    double reaction_plane = 0;    // Direction of flow
    double v2 = 0.05;    // v2 and v3 of distribution
    double v3 = 0;

    int multiplicity = 300; // Mean of multiplicity distribution
    int multiplicity_deviation = 25; // Std dev of multiplicity distribution

    double temperature = 0.5; // GeV, temperature of fireball
    double max_track_pt = 10;

    double jet_reso = 0.2;


    // Initialize objects
    TFile *outfile = new TFile(Form("v2_jetfinding_%i.root", job_id), "RECREATE");

    TRandom *rng = new TRandom3();
    rng->SetSeed(job_id); // Selected by roll of fair die

    // Distribution of tracks in eta-phi
    TH2D *track_loc_dist = new TH2D("track_loc", "Track Location", 25, -1 * max_eta, max_eta, 25, 0, TMath::TwoPi());
    track_loc_dist->GetXaxis()->SetTitle("#eta");
    track_loc_dist->GetYaxis()->SetTitle("#phi");
    track_loc_dist->GetZaxis()->SetTitle("Counts");
    track_loc_dist->SetDirectory(outfile);

    // Distribution of track momentum
    TH1D *track_pt_dist = new TH1D("track_pt", "Track p_{T}", 40, 0, max_track_pt);
    track_pt_dist->GetXaxis()->SetTitle("Track p_{T}");
    track_pt_dist->GetYaxis()->SetTitle("Counts");
    track_pt_dist->SetDirectory(outfile);

    // Distribution of jets in eta-phi
    TH2D *jet_loc_dist = new TH2D("jet_loc", "Jet Location", 25, -1 * max_eta, max_eta, 25, 0, TMath::TwoPi());
    jet_loc_dist->GetXaxis()->SetTitle("#eta");
    jet_loc_dist->GetYaxis()->SetTitle("#phi");
    jet_loc_dist->GetZaxis()->SetTitle("Counts");
    jet_loc_dist->SetDirectory(outfile);

    // Distribution of jet momentum
    TH1D *jet_pt_dist = new TH1D("jet_pt", "Jet p_{T}", 40, 0, 30);
    jet_pt_dist->GetXaxis()->SetTitle("Jet p_{T}");
    jet_pt_dist->GetYaxis()->SetTitle("Counts");
    jet_pt_dist->SetDirectory(outfile);

    // Jet finding 
    std::vector<fastjet::PseudoJet> tracks;
    fastjet::GhostedAreaSpec area_spec = fastjet::GhostedAreaSpec(1);
    fastjet::AreaDefinition jet_area = fastjet::AreaDefinition(fastjet::active_area, area_spec);
    fastjet::JetDefinition jet_definition = fastjet::JetDefinition(fastjet::antikt_algorithm, jet_reso);    

    // Event loop
    for (uint i = 0; i < n; i++) {
        if (i % 1000 == 0) {
            std::cout << "Processed " << i << " events" << std::endl;
        }

        tracks.clear();
        
        // Populate event
        int event_multiplicity = (int)rng->Gaus(multiplicity, multiplicity_deviation);
        for (uint j = 0; j < event_multiplicity; j++) {
            double track_eta = max_eta * (2 * rng->Rndm() - 1);
            double track_phi = rng->Rndm() * TMath::TwoPi();
            double delta_phi = -1 * v2 * sin(2 * track_phi - reaction_plane);   // Add flow
            track_phi += delta_phi;
            if (track_phi < 0) {
                track_phi += TMath::TwoPi();
            }
            if (track_phi > TMath::TwoPi()) {
                track_phi -= TMath::TwoPi();
            }
            track_loc_dist->Fill(track_eta, track_phi);
            
            double track_pt = boltzmann_distribution(rng, temperature, max_track_pt);
            track_pt_dist->Fill(track_pt);
            // Turn into pseudojet and add to collection
            fastjet::PseudoJet track;
            track.reset_PtYPhiM(track_pt, track_eta, track_phi);
            tracks.push_back(track);
        }
        
        // Jet finding
        auto cluster_sequence = fastjet::ClusterSequenceArea(tracks, jet_definition, jet_area);
        std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cluster_sequence.inclusive_jets(10));

        // Populate jet statistics
        for (fastjet::PseudoJet jet : jets) {
            if (abs(jet.eta()) > max_eta - jet_reso) {
                continue;
            } 
            jet_loc_dist->Fill(jet.eta(), jet.phi());
            jet_pt_dist->Fill(jet.pt());
        }
    }
    track_loc_dist->Write();
    track_pt_dist->Write();
    jet_loc_dist->Write();
    jet_pt_dist->Write();

    outfile->Close();
}


// Plot the results of the simulation
void plot() {
    TFile *infile = new TFile("v2_jetfinding.root");
    TH2 *track_loc = infile->Get<TH2D>("track_loc");
    TH2 *jet_loc = infile->Get<TH2D>("jet_loc");

    // Get phi component of distribution
    TH1 *track_phi = track_loc->ProjectionY("track_phi");
    TH1 *jet_phi = jet_loc->ProjectionY("jet_phi");
    
    track_phi->Scale(1.0 / track_phi->GetEntries());
    jet_phi->Scale(1.0 / jet_phi->GetEntries());

    TF1 *v2_fit_track = new TF1("v2_fit_track", "[0] * (1 + [1] * 2 * cos(2*x))", 0, TMath::TwoPi()); // Fit distribution on [0, pi]
    TF1 *v2_fit_jet = new TF1("v2_fit_jet", "[0] * (1 + [1] * 2 * cos(2*x))", 0, TMath::TwoPi()); // Fit distribution on [0, pi]
    v2_fit_track->SetParameters(1, 0);
    v2_fit_jet->SetParameters(1, 0);
    v2_fit_track->SetParNames("offset", "v2");
    v2_fit_jet->SetParNames("offset", "v2");
    track_phi->Fit(v2_fit_track);
    double v2_track = v2_fit_track->GetParameter("v2");
    jet_phi->Fit(v2_fit_jet);
    double v2_jet = v2_fit_jet->GetParameter("v2");

    // Plotting!
    TCanvas *c = new TCanvas("", "", 1000, 1000);
    c->SetLeftMargin(0.15);

    THStack *stack = new THStack();
    track_phi->SetLineColor(kRed);
    v2_fit_track->SetLineColor(kRed);
    jet_phi->SetLineColor(kBlue);
    v2_fit_jet->SetLineColor(kBlue);
    
    stack->Add(track_phi);
    stack->Add(jet_phi);


    TLegend *legend = new TLegend(0.2, 0.12, 0.5, 0.22);
    legend->AddEntry(track_phi, Form("Tracks: v_{2}=%.3f", v2_track));
    legend->AddEntry(jet_phi, Form("Jets: v_{2}=%.3f", v2_jet));
    legend->SetBorderSize(0);
    legend->SetLineWidth(15);
    legend->SetTextSize(0.03);
    
    stack->Draw("nostack e");
    stack->GetXaxis()->SetTitle("#phi");
    stack->GetYaxis()->SetTitle("Frequency");

    v2_fit_track->Draw("same");
    v2_fit_jet->Draw("same");

    legend->Draw();

    TLatex *text = new TLatex();
    text->SetTextSize(0.03);
    text->SetTextFont(42);
    text->DrawLatexNDC(0.2, 0.32, "Combinatorial Jet v_{2} Effects");
    text->DrawLatexNDC(0.2, 0.28, "Assumed Particle v_{2} = 0.05");
    text->DrawLatexNDC(0.2, 0.24, "Anti-k_{T} R=0.2, p_{T}^{#hbox{jet}} > 10 GeV");

    c->SaveAs("toy_v2.png");
}

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cout << "Specify `sim` or `plot`" << std::endl;
        return 1;
    }
    if (strcmp(argv[1], "sim") == 0) {
        if (argc < 3) {
            std::cout << "Specify job id" << std::endl;
            return 1;
        }
        simulation(std::atoi(argv[2]));
        return 0;
    }
    if (strcmp(argv[1], "plot") == 0) {
        plot();
        return 0;
    }
}