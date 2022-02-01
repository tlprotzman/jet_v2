// Testing if I can get the jet reader class working
// Tristan Protzman
// Lehigh University 
// tlp220@lehigh.edu
// December 27, 2021

// ROOT Includes
#include "TROOT.h"
#include "TCanvas.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TMath.h"
#include "TPad.h"
#include "TVector3.h"
#include "THStack.h"
#include "TLegend.h"

// Jetreader Includes
#include "jetreader/reader/reader.h"

// Fastjet Includes
#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/tools/BackgroundEstimatorBase.hh"
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"

// EPD Includes
#include "StEpdUtil/StEpdEpFinder.h"
#include "StEpdUtil/StEpdEpInfo.h"
#include "StPicoEvent/StPicoEpdHit.h"

// C/C++ Libraries
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>


int main(int argc, char **argv) {
    // PARAMETERS
    int DEBUG_LEVEL = 0;
    if (getenv("JR_DEBUG_LEVEL") != nullptr) {
        DEBUG_LEVEL = atoi(getenv("JR_DEBUG_LEVEL"));
    }
    std::cout << "Running at debug level " << DEBUG_LEVEL << std::endl;
    
    // std::string picos_to_read = "test.list"; 
    std::string picos_to_read = "smalltest.list"; 
    // std::string picos_to_read = "/data/star/production_isobar_2018/ReversedFullField/P20ic/2018/083/19083049/st_physics_19083049_raw_1000011.picoDst.root";
    std::string bad_run_list = "";


    // INITIALIZATION
    jetreader::Reader *reader = new jetreader::Reader(picos_to_read);
    // CUTS
    // Are we using a bad run list?
    if (bad_run_list != "") {
        reader->eventSelector()->addBadRuns(bad_run_list);
    }
    
    // Vertex Position
    reader->eventSelector()->setVzRange(-50, 50);
    reader->eventSelector()->setVrMax(0.5);

    // Track cuts
    reader->trackSelector()->setPtMax(30.0);
    reader->trackSelector()->rejectEventOnPtFailure(true);
    // reader->trackSelector()->setDcaMax(trackDCACut);
    // reader->trackSelector()->setNHitsMin(trackNhitCut);
    // reader->trackSelector()->setNHitsFracMin(trackNhitFracCut);

    // Tower cuts
    reader->useHadronicCorrection(true, 1.0);
    reader->useApproximateTrackTowerMatching(true);
    std::string badtowerlist = "";
    // if(removebadtowers){
    //   badtowerlist = "test.csv";
    //   reader->towerSelector()->addBadTowers(badtowerlist);
    // }
    reader->towerSelector()->setEtMax(30.0);
    reader->towerSelector()->rejectEventOnEtFailure(true);

    
    // Initialize Reader
    reader->Init();
    if (DEBUG_LEVEL > 0) {
        std::cout << reader->tree()->GetEntries() << " events considered" << std::endl;
    }    

    
    // Set up jet finding
    double jet_r = 0.3;
    
    double ghost_maxrap = 1;
    fastjet::GhostedAreaSpec area_spec(ghost_maxrap);
    fastjet::AreaDefinition jet_area(fastjet::active_area, area_spec);
    fastjet::AreaDefinition jet_area_background(fastjet::active_area_explicit_ghosts, fastjet::GhostedAreaSpec(ghost_maxrap));

    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, jet_r);
    fastjet::JetDefinition jet_subtraction_jet_def(fastjet::kt_algorithm, jet_r);

    fastjet::Selector jet_selector = fastjet::SelectorAbsRapMax(1) * (!fastjet::SelectorNHardest(2));
    fastjet::JetMedianBackgroundEstimator jet_background_estimator(jet_selector, jet_subtraction_jet_def, jet_area_background);
    fastjet::GridMedianBackgroundEstimator grid_background_estimator(1, 0.5);

    fastjet::Subtractor jet_backgorund_subtractor(&jet_background_estimator);
    fastjet::Subtractor grid_backgorund_subtractor(&grid_background_estimator);

    jet_backgorund_subtractor.set_use_rho_m(true);
    jet_backgorund_subtractor.set_safe_mass(true);

    grid_backgorund_subtractor.set_use_rho_m(true);
    grid_backgorund_subtractor.set_safe_mass(true);

    // Histograms
    TH1D *track_constituent_momentum = new TH1D("track_constituent_momentum", "Track Momentum", 50, 0, 40);
    TH1D *jet_momentum = new TH1D("jet_momentum", "Jet Momentum", 50, 0, 60);
    TH1D *jet_momentum_jet_median_subtracted = new TH1D("jet_momentum_median_subtracted", "Jet Momentum, Median Subtracted", 50, -15, 60);
    TH1D *jet_momentum_grid_median_subtracted = new TH1D("jet_momentum_grid_subtracted", "Jet Momentum, Grid Subtracted", 50, -15, 60);
    std::vector<TH1D*> ep(6); // e_uncorrected, w_uncorrected, e_phi, w_phi, e_phi_psi, w_phi_psi
    const char *ep_hist_names[6] = {"east_uncorrected", 
                                    "west_uncorrected",
                                    "east_phi_corrected",
                                    "west_phi_corrected",
                                    "east_phi_psi_corrected",
                                    "west_phi_psi_corrected",};
    for (uint32_t i = 0; i < 6; i++) {
        ep[i] = new TH1D(ep_hist_names[i], ep_hist_names[i], 30, 0, TMath::TwoPi());
    }
    TH1D *hardcore_jets_pt = new TH1D("hardcore_jets", "Hardcore Jet Pt", 50, -15, 60);
    TH1D *matched_jets = new TH1D("matched_jets", "Matched Jet Pt", 50, -15, 60);
    TH1D *jet_z = new TH1D("jet_z", "Jet Z", 25, 0, 1);
    

    // Set up event plane finding
    StEpdEpFinder *ep_finder = new StEpdEpFinder(1, "StEpdEpFinderCorrectionHistograms_OUTPUT.root", "StEpdEpFinderCorrectionHistograms_INPUT.root");
    ep_finder->SetnMipThreshold(0.3);
    ep_finder->SetMaxTileWeight(3);
    ep_finder->SetEpdHitFormat(2); // for StPicoDst


    // Event Loop
    int processed_events = 0;
    while (reader->next()) {
        processed_events++;
        if (processed_events % 1000 == 0) {
            std::cout << "Processed " << processed_events << " events" << std::endl;
        }

        // Find event plane
        TVector3 primary_vertex = reader->picoDst()->event()->primaryVertex();
        TClonesArray *epd_hits = reader->picoDst()->picoArray(8);
        StEpdEpInfo ep_info = ep_finder->Results(epd_hits, primary_vertex, 0);
        ep[0]->Fill(ep_info.EastRawPsi(1));
        ep[1]->Fill(ep_info.WestRawPsi(1));
        ep[2]->Fill(ep_info.EastPhiWeightedPsi(1));
        ep[3]->Fill(ep_info.WestPhiWeightedPsi(1));
        ep[4]->Fill(ep_info.EastPhiWeightedAndShiftedPsi(1));
        ep[5]->Fill(ep_info.WestPhiWeightedAndShiftedPsi(1));
        


        std::vector<fastjet::PseudoJet> tracks = reader->pseudojets();
        std::vector<fastjet::PseudoJet> hardcore_tracks;
        for (std::vector<fastjet::PseudoJet>::iterator track = tracks.begin(); track != tracks.end(); track++) {    // Let's just play with tpc tracks for now
            jetreader::VectorInfo track_info = track->user_info<jetreader::VectorInfo>();
            if ((!track_info.isPrimary())) {
                tracks.erase(track);
                track--;
            }
            else {
                double track_pt = track->pt();
                track_constituent_momentum->Fill(track_pt);
                if (track_pt > 2) {
                    hardcore_tracks.push_back(*track);
                }
            }
        }

        fastjet::ClusterSequenceArea hc_cs = fastjet::ClusterSequenceArea(hardcore_tracks, jet_def, jet_area);
        fastjet::ClusterSequenceArea cs = fastjet::ClusterSequenceArea(tracks, jet_def, jet_area);

        std::vector<fastjet::PseudoJet> hardcore_jets = fastjet::sorted_by_pt(hc_cs.inclusive_jets());
        std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets()); // TODO inclusive vs exclusive jets

        jet_background_estimator.set_particles(tracks);
        grid_background_estimator.set_particles(tracks);

        for (fastjet::PseudoJet hc_jet : hardcore_jets) {
            hardcore_jets_pt->Fill(hc_jet.pt());
            for (fastjet::PseudoJet jet : jets) {
                if (hc_jet.squared_distance(jet) < 0.3){
                    double max_pt = 0;
                    for (auto constituent : jet.constituents()) {
                        if (constituent.pt() > max_pt) {
                            max_pt = constituent.pt();
                        }
                    }
                    jet_z->Fill(max_pt / jet.pt());
                    if (jet.pt() > 20 && max_pt / jet.pt() < 0.85){
                        jet_momentum->Fill(jet.pt());
                        jet_momentum_jet_median_subtracted->Fill(jet.pt() - jet_background_estimator.rho() * jet.area_4vector().pt());
                        jet_momentum_grid_median_subtracted->Fill(jet.pt() - grid_background_estimator.rho() * jet.area_4vector().pt());
                    }
                }
            }
        }

    }
    std::cout << "Count: " << processed_events << std::endl;
    ep_finder->Finish();

    // Quick plotting
    TCanvas *canvas = new TCanvas("", "", 1000, 1000);
    jet_momentum->Draw("hist");
    gPad->SetLogy();
    canvas->Print("plots/jet_momentum.png");
    jet_momentum_jet_median_subtracted->Draw("hist");
    gPad->SetLogy();
    canvas->Print("plots/jet_momentum_jet_subtracted.png");
    jet_momentum_grid_median_subtracted->Draw("hist");
    gPad->SetLogy();
    canvas->Print("plots/jet_momentum_grid_subtracted.png");
    track_constituent_momentum->Draw("hist");
    gPad->SetLogy();
    canvas->Print("plots/track_constituent_momentum.png");
    hardcore_jets_pt->Draw("hist");
    gPad->SetLogy();
    canvas->Print("plots/hardcore_jet_pt.png");
    jet_z->Draw("hist");
    gPad->SetLogy();
    canvas->Print("plots/jet_z.png");
    TCanvas *ep_canvas = new TCanvas("ep", "", 1000, 1000);
    for (uint32_t i = 0; i < 6; i+=2) {
        THStack *ep_stack = new THStack();
        TLegend *ep_legend = new TLegend();
        
        ep[i]->SetLineColor(kRed);
        ep_stack->Add(ep[i]);
        ep_legend->AddEntry(ep[i], "East");
        
        ep[i+1]->SetLineColor(kBlue);
        ep_stack->Add(ep[i+1]);
        ep_legend->AddEntry(ep[i+1], "West");
        
        ep_stack->Draw("nostack");
        ep_stack->GetXaxis()->SetTitle("Phi");
        ep_stack->GetYaxis()->SetTitle("Counts");
        ep_stack->SetTitle(ep_hist_names[i] + 5);
        ep_legend->Draw();
        ep_canvas->Print(Form("plots/%s.png", ep_hist_names[i] + 5));
        delete ep_stack;
        delete ep_legend;
    }
}

