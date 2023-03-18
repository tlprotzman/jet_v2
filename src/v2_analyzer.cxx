#include "v2_analyzer.h"

#include <TROOT.h>
#include <TChain.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2.h>
#include <TH2D.h>
#include <TH3.h>
#include <TH3D.h>
#include <TFile.h>
#include <TMath.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TLatex.h>
#include <TRandom.h>
#include <TRandom3.h>

#include <time.h>

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <set>

// CONSTRUCTOR
v2_analyzer::v2_analyzer(double jet_resolution, std::string file_list) {
    this->jet_reso = jet_resolution;
    this->file_list = file_list;
    this->load_files();
    this->register_histograms();
    this->event_tree = new Event_Tree(this->chain, "event");
    this->event_tree->readable_tree();
    this->jets = new Jet_Tree(this->chain, "all");
    this->jets->readable_tree();
    this->hc_jets = new Jet_Tree(this->chain, "hc");
    this->hc_jets->readable_tree();

    this->rng = new TRandom3(time(nullptr));    // Seed RNG with time, fix if running on clusters
    this->used_hc_jets = new std::set<int>;

    this->num_events = 0;
    this->event_plane_resolution = 0;
    std::cout << "V2 Analyzer Initialized" << std::endl;
}

// PUBLIC
// Creates histogram map and populates
void v2_analyzer::run_analysis(int n) {
    this->event_loop(n);
    this->calculate_v2();
}

// Saves all the histograms in the map to the specified file
void v2_analyzer::save_histograms(std::string file_name) {
    TFile *outfile = new TFile(file_name.c_str(), "RECREATE");
    for (TH1 *hist : *this->histograms) {
        hist->SetDirectory(outfile);
        hist->Write();
        std::cout << "Saved " << hist->GetName() << std::endl;
    }
    outfile->GetList()->Add(this->jet_v2_no_modulation);  // This is frustrating
    outfile->GetList()->Add(this->jet_v2_extra_modulation);  // This is frustrating
    outfile->GetList()->Add(this->jet_v2_hc);  
    outfile->GetList()->Add(this->jet_v2_all);  
    outfile->GetList()->Add(this->jet_v2_hc_only);  
    this->jet_v2_no_modulation->Write();
    this->jet_v2_extra_modulation->Write();
    this->jet_v2_hc->Write();
    this->jet_v2_all->Write();
    this->jet_v2_hc_only->Write();
    outfile->Close();
    delete outfile;
}

// PRIVATE
// Load trees to analyze and disables all branches
void v2_analyzer::load_files() {
    this->chain = new TChain("jet_data");
    std::ifstream files(this->file_list);
    std::string filePath;
    while (files >> filePath) {
        this->chain->AddFile(filePath.c_str());
    }
    files.close();
    this->chain->SetBranchStatus("*", 0);
    this->chain->SetBranchStatus("*_jet_eta", 1);
    this->chain->SetBranchStatus("*_jet_phi", 1);
    this->chain->SetBranchStatus("*_num_jets", 1);
    this->chain->SetBranchStatus("*_jet_pt", 1);
    this->chain->SetBranchStatus("all_jet_area_pt", 1);
    this->chain->SetBranchStatus("event_ep*", 1);
    this->chain->SetBranchStatus("event_centrality", 1);
}

// Initializes the histogram manager
void v2_analyzer::register_histograms() {
    this->histograms = new std::vector<TH1*>;

    // Jet vs EP hist;
    int n_phi_bins = 4;
    double phi_low = 0;
    double phi_high = TMath::PiOver2();
    double phi_bins[n_phi_bins + 1];
    for (uint32_t i = 0; i < n_phi_bins + 1; i++) {
        phi_bins[i] = phi_low + i * ((phi_high - phi_low) / n_phi_bins);
    }


    this->jet_relative_angle_hc = new TH2D("relative_hc", "relative hc",
                      this->n_pt_bins_v2, this->pt_bins,
                      n_phi_bins, phi_bins);

    this->jet_relative_angle_hc->GetXaxis()->SetTitle("Momentum");
    this->jet_relative_angle_hc->GetYaxis()->SetTitle("Phi");
    this->jet_relative_angle_hc->SetTitle("Jet Azimuthal Angle Relative to Event Plane");
    histograms->push_back(this->jet_relative_angle_hc);

    this->jet_relative_angle_no_background_mod = new TH2D("relative_no_mod", "relative no mod",
                      this->n_pt_bins_v2, this->pt_bins,
                      n_phi_bins, phi_bins);

    this->jet_relative_angle_no_background_mod->GetXaxis()->SetTitle("Momentum");
    this->jet_relative_angle_no_background_mod->GetYaxis()->SetTitle("Phi");
    this->jet_relative_angle_no_background_mod->SetTitle("Jet Azimuthal Angle Relative to Event Plane");
    histograms->push_back(this->jet_relative_angle_no_background_mod);

    this->jet_relative_angle_extra_background_mod = new TH2D("relative_extra_mod", "relative avg extra mod",
                      this->n_pt_bins_v2, this->pt_bins,
                      n_phi_bins, phi_bins);

    this->jet_relative_angle_extra_background_mod->GetXaxis()->SetTitle("Momentum");
    this->jet_relative_angle_extra_background_mod->GetYaxis()->SetTitle("Phi");
    this->jet_relative_angle_extra_background_mod->SetTitle("Jet Azimuthal Angle Relative to Event Plane");
    histograms->push_back(this->jet_relative_angle_extra_background_mod);

    this->jet_relative_angle_all = new TH2D("relative_all", "relative all",
                      this->n_pt_bins_v2, this->pt_bins,
                      n_phi_bins, phi_bins);

    this->jet_relative_angle_all->GetXaxis()->SetTitle("Momentum");
    this->jet_relative_angle_all->GetYaxis()->SetTitle("Phi");
    this->jet_relative_angle_all->SetTitle("Jet Azimuthal Angle Relative to Event Plane");
    histograms->push_back(this->jet_relative_angle_all);

    this->jet_relative_angle_hc_only = new TH2D("relative_hc_only", "relative hc_only",
                      this->n_pt_bins_v2, this->pt_bins,
                      n_phi_bins, phi_bins);

    this->jet_relative_angle_hc_only->GetXaxis()->SetTitle("Momentum");
    this->jet_relative_angle_hc_only->GetYaxis()->SetTitle("Phi");
    this->jet_relative_angle_hc_only->SetTitle("Jet Azimuthal Angle Relative to Event Plane");
    histograms->push_back(this->jet_relative_angle_hc_only);

    // Jet spectra
    int n_pt_bins = 50;
    double pt_low = 0;
    double pt_high = 50;
    this->hc_pt_spectra = new TH1D("hc_pt_spectra", "Hard Core p_{T} Spectra", n_pt_bins, pt_low, pt_high);
    this->hc_pt_spectra->GetXaxis()->SetTitle("p_{T}");
    this->hc_pt_spectra->GetYaxis()->SetTitle("Counts");
    this->histograms->push_back(this->hc_pt_spectra);

    n_pt_bins = 61;
    pt_low = -10;
    this->jet_pt_spectra = new TH1D("jet_pt_spectra", "Jet p_{T} Spectra", n_pt_bins, pt_low, pt_high);
    this->jet_pt_spectra->GetXaxis()->SetTitle("p_{T}");
    this->jet_pt_spectra->GetYaxis()->SetTitle("Counts");
    this->histograms->push_back(this->jet_pt_spectra);

    this->jet_pt_spectra_not_matched = new TH1D("jet_pt_spectra_nomatch", "Jet p_{T} Spectra, no matching requirement", n_pt_bins, pt_low, pt_high);
    this->jet_pt_spectra_not_matched->GetXaxis()->SetTitle("p_{T}");
    this->jet_pt_spectra_not_matched->GetYaxis()->SetTitle("Counts");
    this->histograms->push_back(this->jet_pt_spectra_not_matched);

    // UE Subtraction
    this->ue_subtraction = new TH2D("ue_subtraction", "Underlying Event Subtracted", 25, 0, TMath::PiOver2(), 50, 0, 7);
    this->ue_subtraction->GetXaxis()->SetTitle("#Delta#phi");
    this->ue_subtraction->GetYaxis()->SetTitle("p_{T}");
    this->histograms->push_back(this->ue_subtraction);

    this->subtraction_relative_ep = new TH2D("matchjet_diff_pt_rel_ep", "p_{T} difference relative to EP", 25, 0, TMath::Pi(), 50, -20, 15);
    this->subtraction_relative_ep->GetXaxis()->SetTitle("#Delta #phi");
    this->subtraction_relative_ep->GetYaxis()->SetTitle("p_{T}");
    this->histograms->push_back(this->subtraction_relative_ep);

    this->rho_distribution = new TH1D("rho_dist", "#rho", 50, 0, 25);
    this->rho_distribution->GetXaxis()->SetTitle("#rho (GeV/c / #Omega)");
    this->rho_distribution->GetYaxis()->SetTitle("Counts");
    this->histograms->push_back(this->rho_distribution);

    // Matching performance
    this->matching_performance = new TH3D("matching_performance", "Matching Performance", 50, -1.2 * this->jet_reso, 1.2 * this->jet_reso, 
                                                                                          50, -1.2 * this->jet_reso, 1.2 * this->jet_reso, 
                                                                                          50, -20, 15);
    this->matching_performance->GetXaxis()->SetTitle("#Delta #eta");
    this->matching_performance->GetYaxis()->SetTitle("#Delta #phi");
    this->matching_performance->GetZaxis()->SetTitle("#Delta p_{T}");
    this->histograms->push_back(this->matching_performance);

    // Event Plane Stuff
    this->ep_correlation = new TH2D("ep_correlation", "#Psi_{2} Correlation", 25, 0, TMath::Pi(), 25, 0, TMath::Pi());
    this->ep_correlation->GetXaxis()->SetTitle("East #Psi_{2}");
    this->ep_correlation->GetYaxis()->SetTitle("West #Psi_{2}");
    this->histograms->push_back(this->ep_correlation);

    this->ep_diff = new TH1D("ep_diff", "Event Plane Difference", 50, 0, TMath::Pi());
    this->ep_diff->GetXaxis()->SetTitle("#Psi_{2, east} - #Psi_{2, west}");
    this->histograms->push_back(this->ep_diff);

    this->ep_resolution_cent_binned = new TH2D("ep_reso", "Event Plane Resolution", 16, -0.5, 15.5, 100, 0, TMath::Pi());
    this->ep_resolution_cent_binned->GetXaxis()->SetTitle("Centrality");
    this->ep_resolution_cent_binned->GetYaxis()->SetTitle("#Psi_{2, east} - #Psi_{2, west}");
    this->histograms->push_back(this->ep_resolution_cent_binned);
}

// Determins if a jet matches to a hard core
int v2_analyzer::hardcore_matched(int index) {
    double jet_phi = this->jets->jet_phi[index];
    double jet_eta = this->jets->jet_eta[index];
    int matched = -1;
    double closest = this->jet_reso * this->jet_reso * 1.1; // A value larger that we're searching in
    for (int i = 0; i < this->hc_jets->num_jets; i++) {
        if (this->used_hc_jets->count(i)) {
            // std::cout << "already used" << std::endl;
            continue;
        }
        // Skip hard cores with too little momentum
        if (this->hc_jets->jet_pt[i] < this->hc_pt_min) {
            continue;
        }
        // Fiducal cut
        if (abs(this->hc_jets->jet_eta[i]) > 1 - this->jet_reso) {
            continue;
        }
        double d_phi = abs(jet_phi - this->hc_jets->jet_phi[i]);
        double d_eta = abs(jet_eta - this->hc_jets->jet_eta[i]);
        double dr = d_phi * d_phi + d_eta *d_eta;
        // If dr < jet resolution, we say it is matched
        // std::cout << dr << ":" << closest << std::endl;
        if (dr < closest) {
            // std::cout << "found closer jet" << std::endl;
            matched = i;
            closest = dr;
        }
    }
    // if (matched >= 0) {
    //     std::cout << matched << std::endl;
    // }
    this->used_hc_jets->insert(matched);
    return matched;
}

// Loops over all saved events
void v2_analyzer::event_loop(int n) {
    for (Long64_t i = 0; i < this->chain->GetEntries(); i++) {
        this->chain->GetEvent(i);
        this->used_hc_jets->clear();    // Reset set of used hc jets
        bool match_found = false;
        if (n != 0 && i > n) {
            break;
        }
        if (i % 100000 == 0) {
            printf("Processed %.1fM events\n", i / 1000000.);
        }
        
        // DON'T USE THIS, FASTJET MANAUL SAYS 0 RHO IS NOT UNEXPECTED
        // LEAVING IN SO I DON'T TRY TO ADD IT LATER
        // Get rid of events with no background, something probably went wrong
        // if (this->jets->rho < 0.05) {
        //     continue;
        // }

        double ep_diff_angle = this->event_tree->ep_east - this->event_tree->ep_west;
        if (ep_diff_angle < 0) {
            ep_diff_angle += TMath::TwoPi();
        }
        if (ep_diff_angle > TMath::Pi()) {
            ep_diff_angle -= TMath::Pi();
        }
        this->ep_resolution_cent_binned->Fill(this->event_tree->centrality, ep_diff_angle);


        // Final QA pass
        // Starting with centrality 
        if (this->event_tree->centrality > this->centrality_upper_bound || this->event_tree->centrality < this->centrality_lower_bound) {
            continue;
        }

        this->rho_distribution->Fill(this->jets->rho);

        // Loop over hard cores
        for (int j = 0; j < this->hc_jets->num_jets; j++) {
            // Fiducal cut
            if (abs(this->hc_jets->jet_eta[j]) > 1 - this->jet_reso) {
                continue;
            }
            this->hc_pt_spectra->Fill(this->hc_jets->jet_pt[j]);
            double relative_angle_east = this->hc_jets->jet_phi[j] - this->event_tree->ep_east;
            if (relative_angle_east < 0) {
                relative_angle_east += TMath::TwoPi();
            }
            if (relative_angle_east > TMath::Pi()) {
                relative_angle_east -= TMath::Pi();
            }
            if (relative_angle_east > TMath::PiOver2()) {
                relative_angle_east = TMath::Pi() - relative_angle_east;
            }
            this->jet_relative_angle_hc_only->Fill(this->hc_jets->jet_pt[j], relative_angle_east);
        }
        
        // Looping over jets
        for (int j = 0; j < this->jets->num_jets; j++) {
            // Fiducal cut
            if (abs(this->jets->jet_eta[j]) > 1 - this->jet_reso) {
                continue;
            }
            // Get rid of events with less than 10% the target area
            if (this->jets->jet_area_pt[j]  < 0.1 * TMath::Pi() * this->jet_reso * this->jet_reso) {
                continue;
            }
            // Calculate \Delta\Phi
            double avg_angle = atan2 ( 0.5 * (sin(this->event_tree->ep_east) + sin(this->event_tree->ep_west)), 0.5 * (cos(this->event_tree->ep_east) + cos(this->event_tree->ep_west)));
            double relative_angle_east = this->jets->jet_phi[j] - this->event_tree->ep_east;
            if (relative_angle_east < 0) {
                relative_angle_east += TMath::TwoPi();
            }
            if (relative_angle_east > TMath::Pi()) {
                relative_angle_east -= TMath::Pi();
            }
            if (relative_angle_east > TMath::PiOver2()) {
                relative_angle_east = TMath::Pi() - relative_angle_east;
            }
            double relative_angle_avg = this->jets->jet_phi[j] - avg_angle;
            if (relative_angle_avg < 0) {
                relative_angle_avg += TMath::TwoPi();
            }
            if (relative_angle_avg > TMath::Pi()) {
                relative_angle_avg -= TMath::Pi();
            }
            if (relative_angle_east > TMath::PiOver2()) {    // Fold into {0, pi/2}
                relative_angle_east = TMath::Pi() - relative_angle_east;
            }

            // Calculate and subtract underlying event
            double background = this->jets->rho * (1 +  this->assumed_v2 * cos(2 * relative_angle_east)) * this->jets->jet_area_pt[j];
            double background_no_mod = this->jets->rho * (1 +  0.02 * cos(2 * relative_angle_east)) * this->jets->jet_area_pt[j];
            double background_extra_mod = this->jets->rho * (1 +  0.06 * cos(2 * relative_angle_east)) * this->jets->jet_area_pt[j];
            this->jet_pt_spectra_not_matched->Fill(this->jets->jet_pt[j] - background);
            this->jet_relative_angle_all->Fill(this->jets->jet_pt[j] - background, relative_angle_east);

            
            // Require hardcore matching
            int matched_index = hardcore_matched(j);
            if (matched_index == -1) {
                continue;
            }
            match_found = true;
            // this->avg_cos
            this->jet_pt_spectra->Fill(this->jets->jet_pt[j] - background);
            this->jet_relative_angle_hc->Fill(this->jets->jet_pt[j] - background, relative_angle_east);
            this->jet_relative_angle_no_background_mod->Fill(this->jets->jet_pt[j] - background_no_mod, relative_angle_east);
            this->jet_relative_angle_extra_background_mod->Fill(this->jets->jet_pt[j] - background_extra_mod, relative_angle_east);
            this->ue_subtraction->Fill(relative_angle_east, background);
            this->matching_performance->Fill(this->jets->jet_eta[j] - this->hc_jets->jet_eta[matched_index],
                                             this->jets->jet_phi[j] - this->hc_jets->jet_phi[matched_index],
                                             (this->jets->jet_pt[j]) - this->hc_jets->jet_pt[matched_index]);
            this->subtraction_relative_ep->Fill(relative_angle_east, (this->jets->jet_pt[j] - background) - this->hc_jets->jet_pt[matched_index]);

        }
        // Event-wise quantities
        if (match_found) {
            this->num_events++;
            this->ep_correlation->Fill(this->event_tree->ep_east, this->event_tree->ep_west);
            this->event_plane_resolution += cos(2 * ep_diff_angle);
            this->ep_diff->Fill(ep_diff_angle);
        }
        // Calculate final ep resolution
    }
    this->event_plane_resolution /= this->num_events;
    this->event_plane_resolution = sqrt(this->event_plane_resolution) * sqrt(2);
    std::cout << "Number of Events: " << this->num_events << std::endl;
    std::cout << "Event Plane Resolution: " << this->event_plane_resolution << std::endl;
}

void v2_analyzer::calculate_v2() {
    this->jet_v2_hc = new TGraphErrors(this->n_pt_bins_v2);
    this->jet_v2_hc->SetName("v2_hc");
    this->jet_v2_hc->SetEditable(false);   // sigh...
    this->jet_v2_hc->GetXaxis()->SetTitle("p_{T}");
    this->jet_v2_hc->GetYaxis()->SetTitle("v_{2}^{ch jet}");

    for (int i = 1; i <= this->jet_relative_angle_hc->GetNbinsX(); i++) {
        this->jet_relative_angle_hc->GetXaxis()->SetRange(i, i);
        double pt = this->jet_relative_angle_hc->GetXaxis()->GetBinCenter(i);
        double pt_error = this->jet_relative_angle_hc->GetXaxis()->GetBinWidth(i) / 2;
        TH1 *dphi = this->jet_relative_angle_hc->ProjectionY(Form("relative_angle_hc_%d", i));
        dphi->SetTitle(Form("Jet Yield, %.1f < p_{T}^{reco} < %.1f", pt - pt_error, pt + pt_error));
        TF1 *v2_fit = new TF1(Form("v2_bin_%d", i), "[0] * (1 + [1] * 2 * cos(2*x))", 0, TMath::TwoPi()); // Fit distribution on [0, pi]
        v2_fit->SetParameters(1, 0);
        v2_fit->SetParNames("offset", "v2");
        dphi->Fit(v2_fit, "", "", 0, TMath::PiOver2());
        this->histograms->push_back(dphi);
        double v2 = v2_fit->GetParameter("v2");
        double v2_error = v2_fit->GetParError(v2_fit->GetParNumber("v2"));
        v2 /= this->event_plane_resolution;
        v2_error /= this->event_plane_resolution;
        

        this->jet_v2_hc->SetPoint(i - 1, pt, v2);
        this->jet_v2_hc->SetPointError(i - 1, pt_error, v2_error);
    }

    this->jet_v2_all = new TGraphErrors(this->n_pt_bins_v2);
    this->jet_v2_all->SetName("v2_all");
    this->jet_v2_all->SetEditable(false);   // sigh...
    this->jet_v2_all->GetXaxis()->SetTitle("p_{T}");
    this->jet_v2_all->GetYaxis()->SetTitle("v_{2}^{ch jet}");

    for (int i = 1; i <= this->jet_relative_angle_all->GetNbinsX(); i++) {
        this->jet_relative_angle_all->GetXaxis()->SetRange(i, i);
        double pt = this->jet_relative_angle_all->GetXaxis()->GetBinCenter(i);
        double pt_error = this->jet_relative_angle_all->GetXaxis()->GetBinWidth(i) / 2;
        TH1 *dphi = this->jet_relative_angle_all->ProjectionY(Form("relative_angle_all_%d", i));
        dphi->SetTitle(Form("Jet Yield, %.1f < p_{T}^{reco} < %.1f", pt - pt_error, pt + pt_error));
        TF1 *v2_fit = new TF1(Form("v2_bin_%d", i), "[0] * (1 + [1] * 2 * cos(2*x))", 0, TMath::TwoPi()); // Fit distribution on [0, pi]
        v2_fit->SetParameters(1, 0);
        v2_fit->SetParNames("offset", "v2");
        dphi->Fit(v2_fit, "", "", 0, TMath::PiOver2());
        this->histograms->push_back(dphi);
        double v2 = v2_fit->GetParameter("v2");
        double v2_error = v2_fit->GetParError(v2_fit->GetParNumber("v2"));
        v2 /= this->event_plane_resolution;
        v2_error /= this->event_plane_resolution;
        

        this->jet_v2_all->SetPoint(i - 1, pt, v2);
        this->jet_v2_all->SetPointError(i - 1, pt_error, v2_error);
    }

    this->jet_v2_hc_only = new TGraphErrors(this->n_pt_bins_v2);
    this->jet_v2_hc_only->SetName("v2_hc_only");
    this->jet_v2_hc_only->SetEditable(false);   // sigh...
    this->jet_v2_hc_only->GetXaxis()->SetTitle("p_{T}");
    this->jet_v2_hc_only->GetYaxis()->SetTitle("v_{2}^{ch jet}");

    for (int i = 1; i <= this->jet_relative_angle_hc_only->GetNbinsX(); i++) {
        this->jet_relative_angle_hc_only->GetXaxis()->SetRange(i, i);
        double pt = this->jet_relative_angle_hc_only->GetXaxis()->GetBinCenter(i);
        double pt_error = this->jet_relative_angle_hc_only->GetXaxis()->GetBinWidth(i) / 2;
        TH1 *dphi = this->jet_relative_angle_hc_only->ProjectionY(Form("relative_angle_hc_only_%d", i));
        dphi->SetTitle(Form("Jet Yield, %.1f < p_{T}^{reco} < %.1f", pt - pt_error, pt + pt_error));
        TF1 *v2_fit = new TF1(Form("v2_bin_%d", i), "[0] * (1 + [1] * 2 * cos(2*x))", 0, TMath::TwoPi()); // Fit distribution on [0, pi]
        v2_fit->SetParameters(1, 0);
        v2_fit->SetParNames("offset", "v2");
        dphi->Fit(v2_fit, "", "", 0, TMath::PiOver2());
        this->histograms->push_back(dphi);
        double v2 = v2_fit->GetParameter("v2");
        double v2_error = v2_fit->GetParError(v2_fit->GetParNumber("v2"));
        v2 /= this->event_plane_resolution;
        v2_error /= this->event_plane_resolution;
        

        this->jet_v2_hc_only->SetPoint(i - 1, pt, v2);
        this->jet_v2_hc_only->SetPointError(i - 1, pt_error, v2_error);
    }

    this->jet_v2_no_modulation = new TGraphErrors(this->n_pt_bins_v2);
    this->jet_v2_no_modulation->SetName("v2_no_mod");
    this->jet_v2_no_modulation->SetEditable(false);   // sigh...
    this->jet_v2_no_modulation->GetXaxis()->SetTitle("p_{T}");
    this->jet_v2_no_modulation->GetYaxis()->SetTitle("v_{2}^{ch jet}");

    for (int i = 1; i <= this->jet_relative_angle_no_background_mod->GetNbinsX(); i++) {
        this->jet_relative_angle_no_background_mod->GetXaxis()->SetRange(i, i);
        TH1 *dphi = this->jet_relative_angle_no_background_mod->ProjectionY(Form("relative_angle_no_mod_%d", i));
        TF1 *v2_fit = new TF1(Form("v2_bin_%d", i), "[0] * (1 + [1] * 2 * cos(2*x))", 0, TMath::TwoPi()); // Fit distribution on [0, pi]
        v2_fit->SetParameters(1, 0);
        v2_fit->SetParNames("offset", "v2");
        dphi->Fit(v2_fit, "", "", 0, TMath::PiOver2());
        this->histograms->push_back(dphi);
        double v2 = v2_fit->GetParameter("v2");
        double v2_error = v2_fit->GetParError(v2_fit->GetParNumber("v2"));
        v2 /= this->event_plane_resolution;
        v2_error /= this->event_plane_resolution;
        
        double pt = this->jet_relative_angle_no_background_mod->GetXaxis()->GetBinCenter(i);
        double pt_error = this->jet_relative_angle_no_background_mod->GetXaxis()->GetBinWidth(i) / 2;
        TLatex *pt_range = new TLatex(0.2, 0.2, Form("%.0d<p_{T}<%.0d", pt - pt_error, pt + pt_error));
        TLatex *fit_parameters = new TLatex(0.2, 0.1, Form("%d * (1 + %d * 2 * cos(2*#Delta#phi))", v2_fit->GetParameter("offset"), v2));

        this->jet_v2_no_modulation->SetPoint(i - 1, pt, v2);
        this->jet_v2_no_modulation->SetPointError(i - 1, pt_error, v2_error);
    }

    this->jet_v2_extra_modulation = new TGraphErrors(this->n_pt_bins_v2);
    this->jet_v2_extra_modulation->SetName("v2_extra_mod");
    this->jet_v2_extra_modulation->SetEditable(false);   // sigh...
    this->jet_v2_extra_modulation->GetXaxis()->SetTitle("p_{T}");
    this->jet_v2_extra_modulation->GetYaxis()->SetTitle("v_{2}^{ch jet}");

    for (int i = 1; i <= this->jet_relative_angle_extra_background_mod->GetNbinsX(); i++) {
        this->jet_relative_angle_extra_background_mod->GetXaxis()->SetRange(i, i);
        TH1 *dphi = this->jet_relative_angle_extra_background_mod->ProjectionY(Form("relative_angle_extra_mod_%d", i));
        TF1 *v2_fit = new TF1(Form("v2_bin_%d", i), "[0] * (1 + [1] * 2 * cos(2*x))", 0, TMath::TwoPi()); // Fit distribution on [0, pi]
        v2_fit->SetParameters(1, 0);
        v2_fit->SetParNames("offset", "v2");
        dphi->Fit(v2_fit, "", "", 0, TMath::PiOver2());
        this->histograms->push_back(dphi);
        double v2 = v2_fit->GetParameter("v2");
        double v2_error = v2_fit->GetParError(v2_fit->GetParNumber("v2"));
        v2 /= this->event_plane_resolution;
        v2_error /= this->event_plane_resolution;
        
        double pt = this->jet_relative_angle_extra_background_mod->GetXaxis()->GetBinCenter(i);
        double pt_error = this->jet_relative_angle_extra_background_mod->GetXaxis()->GetBinWidth(i) / 2;
        TLatex *pt_range = new TLatex(0.2, 0.2, Form("%.0d<p_{T}<%.0d", pt - pt_error, pt + pt_error));
        TLatex *fit_parameters = new TLatex(0.2, 0.1, Form("%d * (1 + %d * 2 * cos(2*#Delta#phi))", v2_fit->GetParameter("offset"), v2));

        this->jet_v2_extra_modulation->SetPoint(i - 1, pt, v2);
        this->jet_v2_extra_modulation->SetPointError(i - 1, pt_error, v2_error);
    }

    // this->jet_v2_avg_ep = new TGraphErrors(this->n_pt_bins_v2);
    // this->jet_v2_avg_ep->SetName("v2_avg");
    // this->jet_v2_avg_ep->SetEditable(false);   // sigh...
    // this->jet_v2_avg_ep->GetXaxis()->SetTitle("p_{T}");
    // this->jet_v2_avg_ep->GetYaxis()->SetTitle("v_{2}^{ch jet}");

    // for (int i = 1; i <= this->jet_relative_angle_avg_ep->GetNbinsX(); i++) {
    //     this->jet_relative_angle_avg_ep->GetXaxis()->SetRange(i, i);
    //     TH1 *dphi = this->jet_relative_angle_avg_ep->ProjectionY(Form("relative_angle_%d", i));
    //     TF1 *v2_fit = new TF1(Form("v2_bin_%d", i), "[0] * (1 + [1] * 2 * cos(2*x))", 0, TMath::TwoPi()); // Fit distribution on [0, pi]
    //     v2_fit->SetParameters(1, 0);
    //     v2_fit->SetParNames("offset", "v2");
    //     dphi->Fit(v2_fit, "", "", 0, TMath::PiOver2());
    //     this->histograms->push_back(dphi);
    //     double v2 = v2_fit->GetParameter("v2");
    //     double v2_error = v2_fit->GetParError(v2_fit->GetParNumber("v2"));
    //     v2 /= this->event_plane_resolution;
    //     v2_error /= this->event_plane_resolution;
        
    //     double pt = this->jet_relative_angle_avg_ep->GetXaxis()->GetBinCenter(i);
    //     double pt_error = this->jet_relative_angle_avg_ep->GetXaxis()->GetBinWidth(i) / 2;
    //     TLatex *pt_range = new TLatex(0.2, 0.2, Form("%.0d<p_{T}^{reco}<%.0d", pt - pt_error, pt + pt_error));
    //     TLatex *fit_parameters = new TLatex(0.2, 0.1, Form("%d * (1 + %d * 2 * cos(2*#Delta#phi))", v2_fit->GetParameter("offset"), v2));

    //     this->jet_v2_avg_ep->SetPoint(i - 1, pt, v2);
    //     this->jet_v2_avg_ep->SetPointError(i - 1, pt_error, v2_error);
    // }
}

