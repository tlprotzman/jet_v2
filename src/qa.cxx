// Testing if I can get the jet reader class working
// Tristan Protzman
// Lehigh University 
// tlp220@lehigh.edu
// December 27, 2021

// Comment your damn code Tristan!!!

// ROOT Includes
#include <TROOT.h>
#include <TCanvas.h>
#include <TClonesArray.h>
#include <TH1.h>
#include <TMath.h>
#include <TPad.h>
#include <TVector3.h>
#include <THStack.h>
#include <TLegend.h>
#include <TTree.h>

// Jetreader Includes
#include <jetreader/reader/reader.h>

// Fastjet Includes
#include <fastjet/AreaDefinition.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/tools/BackgroundEstimatorBase.hh>
#include <fastjet/tools/GridMedianBackgroundEstimator.hh>
#include <fastjet/tools/JetMedianBackgroundEstimator.hh>
#include <fastjet/tools/Subtractor.hh>

// EPD Includes
#include <StEpdUtil/StEpdEpFinder.h>
#include <StEpdUtil/StEpdEpInfo.h>
#include <StPicoEvent/StPicoEpdHit.h>

// C/C++ Libraries
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>

// My includes
#include "setup.h"
#include "isobar_triggers.h"
#include "event_tree.h"
#include "jet_tree.h"
#include "jet_helper.h"
#include "particle_anisotropy.h"
#include "angle_tools.h"

// double JET_PT_CUT = 5;
double NMIP_MIN = 0.3;
double NMIP_MAX = 3;

// Calculate the sum of nMips in the EPD.  Side 0 is east, 1 is west, 2 is both (Default both).
double epd_mult(TClonesArray *epd_hits, int side=2);
bool pileup_cut(int charged_particles, int tofmatch, int tofmult, int refmult3);

int main(int argc, char **argv) {
    double jet_radius = 0.2;
    // ROOT::EnableImplicitMT();
    // Set up QA manager, controls processes, cuts, io, etc...
    QA_Manager *manager = new QA_Manager(argc, argv);    

    
    // Trees
    TFile *outfile = new TFile(Form("%sout.root", manager->job_id.c_str()), "RECREATE");
    TTree *jet_data = new TTree("jet_data", "Jet Data");
    jet_data->SetDirectory(outfile);
    
    Event_Tree *event_tree = new Event_Tree(jet_data, "event");
    event_tree->writeable_tree();
    Jet_Tree *hardcore_jet_tree_1 = new Jet_Tree(jet_data, "hc_1");
    hardcore_jet_tree_1->writeable_tree();
    Jet_Tree *hardcore_jet_tree_2 = new Jet_Tree(jet_data, "hc_2");
    hardcore_jet_tree_2->writeable_tree();
    Jet_Tree *hardcore_jet_tree_3 = new Jet_Tree(jet_data, "hc_3");
    hardcore_jet_tree_3->writeable_tree();
    Jet_Tree *jet_tree = new Jet_Tree(jet_data, "all");
    jet_tree->writeable_tree();
      
    
    // Set up jet finding
    Jet_Helper *jet_helper = new Jet_Helper(jet_radius);
    Jet_Helper *hardcore_jet_helper = new Jet_Helper(jet_radius);

    // Histograms
    qa_histograms qa_hist;
    ep_histograms ep_hist;
    setup_histograms(&qa_hist, &ep_hist);

    

    // Set up event plane finding
    StEpdEpFinder *ep_finder = new StEpdEpFinder(17, "StEpdEpFinderCorrectionHistograms_OUTPUT.root", "StEpdEpFinderCorrectionHistograms_INPUT.root");
    ep_finder->SetnMipThreshold(NMIP_MIN);
    ep_finder->SetMaxTileWeight(NMIP_MAX);
    ep_finder->SetEpdHitFormat(2); // for StPicoDst

    // Trigger Identification
    isobar_triggers trigger_names; 

    TFile *eventout_file = new TFile("event_display_out.root", "RECREATE");
    

    long possible_tracks = 0;
    long removed_tracks = 0;

    // Event Loop
    int processed_events = 0;
    while (manager->reader->next()) {
        processed_events++;
        if (processed_events % 1000 == 0) {
            std::cout << "Processed " << processed_events << " events" << std::endl;
        }
        if (manager->max_events_set && processed_events >= manager->max_events) {
            break;
        }

        StPicoEvent *event = manager->reader->picoDst()->event();

        // Event Info
        bool is_minbias = false;    // min bias trigger
        bool is_bht1_vpd30 = false; // High tower trigger
        event_tree->num_triggers = 0;
        for (auto t : event->triggerIds()) {
            if (trigger_names.triggers[t] == trigger_names.minbias) {
                is_minbias = true;
            } else if (trigger_names.triggers[t] == trigger_names.bht1_vpd30) {
                is_bht1_vpd30 = true;
            }
            event_tree->triggers[event_tree->num_triggers] = t;
            event_tree->num_triggers++;
            if (event_tree->num_entries <= event_tree->num_triggers) {
                std::cerr << "WARNING: NOT ENOUGH SPACE FOR ALL TRIGGERS" << std::endl;
            }
        }

        // Find event plane
        if (manager->only_ep_finding && !is_minbias) {   // Only use minbias events to generate the EPD corrections
            continue;
        }

        if (manager->reader->centrality16() == -1) { // skip runs we can't determine centrality for
            continue;
        }

        // Pileup cut
        // if (!pileup_cut(event->numberOfPrimaryTracks(), event->nBTOFMatch(), event->btofTrayMultiplicity(), event->refMult3())) {    // This needs to be done properly...
        //     continue;
        // }
        
        event_tree->run_number = event->runId();

        TVector3 primary_vertex = event->primaryVertex();
        event_tree->vx = event->primaryVertex().X();
        event_tree->vy = event->primaryVertex().Y();
        event_tree->vz = event->primaryVertex().Z();
        event_tree->vpd_vz = event->vzVpd();
        event_tree->refmult3 = event->refMult3();
        event_tree->tofmatch = event->nBTOFMatch();
        event_tree->tofmult = event->btofTrayMultiplicity();
        event_tree->bbc_east_rate = event->bbcEastRate();
        event_tree->bbc_west_rate = event->bbcWestRate();
        event_tree->centrality = manager->reader->centrality16();

        qa_hist.vz->Fill(event_tree->vz, event->vzVpd());
        qa_hist.vr->Fill(event_tree->vx, event_tree->vy);
        qa_hist.pileup->Fill(event_tree->refmult3, event_tree->tofmult);
        qa_hist.tofmult->Fill(event_tree->tofmult);
        qa_hist.bbc_rate->Fill(event_tree->bbc_east_rate, event_tree->bbc_west_rate);
        qa_hist.centrality->Fill(event_tree->centrality);
        
       
        TClonesArray *epd_hits = manager->reader->picoDst()->picoArray(8);
        double nMips_sum = epd_mult(epd_hits);
        qa_hist.nMips->Fill(event_tree->refmult3, nMips_sum); // save nmips for jet events too?
        StEpdEpInfo ep_info = ep_finder->Results(epd_hits, primary_vertex, event_tree->centrality);
        ep_hist.east_uncorrected->Fill(ep_info.EastRawPsi(2));
        ep_hist.west_uncorrected->Fill(ep_info.WestRawPsi(2));
        ep_hist.east_phi_corrected->Fill(ep_info.EastPhiWeightedPsi(2));
        ep_hist.west_phi_corrected->Fill(ep_info.WestPhiWeightedPsi(2));
        ep_hist.east_phi_psi_corrected->Fill(ep_info.EastPhiWeightedAndShiftedPsi(2));
        ep_hist.west_phi_psi_corrected->Fill(ep_info.WestPhiWeightedAndShiftedPsi(2));
        ep_hist.ep_correlation->Fill(ep_info.EastPhiWeightedAndShiftedPsi(2), ep_info.WestPhiWeightedAndShiftedPsi(2));
        ep_hist.epd_resolution->Fill(ep_info.EastPhiWeightedAndShiftedPsi(2) - ep_info.WestPhiWeightedAndShiftedPsi(2), manager->reader->centrality16());
        event_tree->ep_east = ep_info.EastPhiWeightedAndShiftedPsi(2);
        event_tree->ep_west = ep_info.WestPhiWeightedAndShiftedPsi(2);
        if (manager->only_ep_finding) {
            // jet_data->Fill();
            continue;
        }
        
        // Calculate charged particle v2 in min-bias data
        if (is_minbias) {
            std::vector<fastjet::PseudoJet> tracks = manager->reader->pseudojets();
            for (std::vector<fastjet::PseudoJet>::iterator track = tracks.begin(); track != tracks.end(); track++) {
                double angle = smallest_angle(track->phi(), event_tree->ep_east);
                qa_hist.particle_v2->Fill(angle, track->pt(), event_tree->centrality);
            }
        }

        // Find jets in events with a high tower trigger
        if (is_bht1_vpd30) {
            std::vector<fastjet::PseudoJet> tracks = manager->reader->pseudojets();
            std::vector<fastjet::PseudoJet> hardcore_tracks_1;
            std::vector<fastjet::PseudoJet> hardcore_tracks_2;
            std::vector<fastjet::PseudoJet> hardcore_tracks_3;
            TH2 *event_display = new TH2D(Form("event_display_%i", processed_events), "Event Display", 25, -1, 1, 25, 0, TMath::TwoPi());
            event_display->GetXaxis()->SetTitle("#eta");
            event_display->GetXaxis()->SetTitle("#phi");
            event_display->GetZaxis()->SetTitle("p_{T}");
            event_display->SetDirectory(eventout_file);
            double max_pt = 0;
            for (std::vector<fastjet::PseudoJet>::iterator track = tracks.begin(); track != tracks.end(); track++) {    // Let's just play with tpc tracks for now
                jetreader::VectorInfo track_info = track->user_info<jetreader::VectorInfo>();
                if (track_info.isBemcTower()) {
                    std::cout << "found tower " << track_info.towerId() << std::endl;
                    
                }
                if ((!track_info.isPrimary())) {
                    tracks.erase(track);
                    track--;
                    continue;
                }
                possible_tracks++;
                if (manager->tracking_efficiency_study && manager->rng->Rndm() < manager->tracking_efficiency_cut) {
                    // std::cerr << "Throwing out track!" << std::endl;
                    removed_tracks++;
                    tracks.erase(track);
                    track--;
                    continue;;
                }
                
                double track_pt = track->pt();
                if (track_pt > max_pt) {
                    max_pt = track_pt;
                }
                qa_hist.track_momentum->Fill(track_pt);
                qa_hist.track_loc->Fill(track->eta(), track->phi());
                event_display->Fill(track->eta(), track->phi(), track->pt());
                if (track_pt > 1) {
                    hardcore_tracks_1.push_back(*track);
                }
                if (track_pt > 2) {
                    hardcore_tracks_2.push_back(*track);
                } 
                if (track_pt > 3) {
                    hardcore_tracks_3.push_back(*track);
                }
            }
            if (max_pt < 8) {
                delete event_display;
            }
            // std::cout << "Regular track count: " << tracks.size() << std::endl;
            // std::cout << "Hardcore track count: " << hardcore_tracks.size() << std::endl;


            std::vector<fastjet::PseudoJet> all_jets = jet_helper->find_jets(tracks);
            jet_helper->set_background_particles(tracks);
            jet_helper->fill_jet_tree(all_jets, jet_tree);
            // jet_helper->fill_jet_tree_particles(tracks, jet_tree);   // Hacky to calculate charged particle v2
            particle_anisotropy flow_description;
            flow_description.rho = -999;
            flow_description.v2 = -999;
            flow_description.v3 = -999;
            // std::cout << "num jets: " << all_jets.size() << std::endl;
            // if (all_jets.size() > 0) {
                // flow_description.calculate_v2(tracks, all_jets[0], (ep_info.EastPhiWeightedAndShiftedPsi(1) + ep_info.WestPhiWeightedAndShiftedPsi(1)) / 2.0, flow_description, manager->reader->centrality16());
            // }
            // printf("Rho: %.4f\tV2: %.4f\tV3: %.04f\n", flow_description.rho, flow_description.v2, flow_description.v3);
            qa_hist.rho->Fill(flow_description.rho, manager->reader->centrality16());
            qa_hist.v2->Fill(flow_description.v2, manager->reader->centrality16());
            qa_hist.v3->Fill(flow_description.v3, manager->reader->centrality16());
            
            hardcore_jet_tree_1->num_jets = 0;
            // std::vector<fastjet::PseudoJet> hardcore_jets_1 = hardcore_jet_helper->find_jets(hardcore_tracks_1);
            // hardcore_jet_helper->set_background_particles(hardcore_tracks_1);
            // hardcore_jet_helper->fill_jet_tree(hardcore_jets_1, hardcore_jet_tree_1);

            std::vector<fastjet::PseudoJet> hardcore_jets_2 = hardcore_jet_helper->find_jets(hardcore_tracks_2);
            hardcore_jet_helper->set_background_particles(hardcore_tracks_2);
            hardcore_jet_helper->fill_jet_tree(hardcore_jets_2, hardcore_jet_tree_2);

            hardcore_jet_tree_3->num_jets = 0;
            // std::vector<fastjet::PseudoJet> hardcore_jets_3 = hardcore_jet_helper->find_jets(hardcore_tracks_3);
            // hardcore_jet_helper->set_background_particles(hardcore_tracks_3);
            // hardcore_jet_helper->fill_jet_tree(hardcore_jets_3, hardcore_jet_tree_3);
            

            // Calculate delta pt
            int tries = 0;
            bool valid = true;
            fastjet::PseudoJet *axis = new fastjet::PseudoJet();
            do {
                double r_eta = manager->rng->Rndm() * 2 - 1;
                double r_phi = manager->rng->Rndm() * TMath::TwoPi();
                axis->reset_momentum_PtYPhiM(1, r_eta, r_phi);
                for (fastjet::PseudoJet jet : hardcore_jets_2) {
                    if (jet.pt() < 5) {
                        continue;
                    }
                    if (axis->delta_R(jet) < jet_radius) {
                        valid = false;
                    }
                }
                tries++;
            } while (valid && tries < 10);
            if (valid) {
                double total_pt = 0;
                for (fastjet::PseudoJet track : tracks) {
                    if (axis->delta_R(track) < jet_radius) {
                        total_pt += track.pt();
                    }
                }
                double pt_diff = total_pt - (TMath::Pi() * jet_radius * jet_radius * jet_tree->rho);
                double delta_ep = axis->phi() - event_tree->ep_east;
                if (delta_ep < 0) {
                    delta_ep += TMath::TwoPi();
                }
                if (jet_tree->rho != 0) {
                    qa_hist.delta_pt->Fill(delta_ep, pt_diff);
                }
            } else {
                std::cerr << "Cannot find valid angle for delta pt calculation" << std::endl;
            }
            
            
            
            // std::cout << "Regular jet cout: " << jet_tree->num_jets << std::endl;
            // std::cout << "Hardcore jet cout: " << hardcore_jet_tree->num_jets << std::endl;

            if (jet_tree->num_jets) {
                event_tree->fill_tree();
            }
        }
    }
    eventout_file->Write();
    eventout_file->Close();

    std::cout << "Possible Tracks: " << possible_tracks << std::endl;
    std::cout << "Removed Tracks: " << removed_tracks << std::endl;
    std::cout << "Removed Fraction: " << (double)removed_tracks / (double)possible_tracks << std::endl;

    std::cout << "Count: " << processed_events << std::endl;
    ep_finder->Finish();


    // Write data to file
    outfile->cd();
    jet_data->Write();
    save_histograms(&qa_hist, &ep_hist, outfile);
    outfile->Close();
    // delete jet_data;
    // cleanup(&datum);
    delete manager;
    delete jet_helper;
    delete hardcore_jet_helper;
}


double epd_mult(TClonesArray *epd_hits, int side) {   // Is this really as slow as it's behaving?
    double nMip_sum = 0;
    for (auto h : *epd_hits) {
        StPicoEpdHit *hit = (StPicoEpdHit*)h;
        if (side == 2 || side == (hit->id() < 0 ? 0 : 1)) {     // 0 is east, 1 is west?  Check this
            double nMips = hit->nMIP();
            if (nMips < NMIP_MIN) {
                nMips = 0;
            }
            if (nMips > NMIP_MAX) {
                nMips = NMIP_MAX;
            }
            nMip_sum += nMips;
        }
    }
    return nMip_sum;
}

// Pileup Cut - stolen shamelessly from https://drupal.star.bnl.gov/STAR/system/files/Isobar_Run18_Step2_QA_Oct14_0_0.pdf
// Returns true if it passes the cut, false otherwise
bool pileup_cut(int charged_particles, int tofmatch, int tofmult, int refmult3) {
    double p0_low = -13.8407;
    double p1_low = 1.00334;
    double p2_low = 0.000421093;
    double p3_low = -1.88309e-06;
    double p4_low = 2.27559e-09;

    double p0_high = 10.4218;
    double p1_high = 1.84005;
    double p2_high = -0.00289939;
    double p3_high = 1.01996e-05;
    double p4_high = -1.472e-08;

    double val_low = p0_low;
    val_low += p1_low * pow(tofmatch, 1);
    val_low += p2_low * pow(tofmatch, 2);
    val_low += p3_low * pow(tofmatch, 3);
    val_low += p4_low * pow(tofmatch, 4);

    double val_high = p0_high;
    val_high += p1_high * pow(tofmatch, 1);
    val_high += p2_high * pow(tofmatch, 2);
    val_high += p3_high * pow(tofmatch, 3);
    val_high += p4_high * pow(tofmatch, 4);

    // std::cout << "Low: " << val_low << "\tHigh: " << val_high << "\n";

    // basic cuts on refmult3 vs tofmult
    bool passA = tofmult > -30 + 1.95 * refmult3 && tofmult < 25 + 2.2 * refmult3;
    // bool passB = charged_particles > val_low && charged_particles < val_high;;
    bool passB = true;
    // std::cout << passA << "\t" << passB << std::endl;
    return passA && passB;
}
