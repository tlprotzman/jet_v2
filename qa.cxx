// Testing if I can get the jet reader class working
// Tristan Protzman
// Lehigh University 
// tlp220@lehigh.edu
// December 27, 2021

// Comment your damn code Tristan!!!

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
#include "TTree.h"

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

// double JET_PT_CUT = 5;
double NMIP_MIN = 0.3;
double NMIP_MAX = 3;

// Calculate the sum of nMips in the EPD.  Side 0 is east, 1 is west, 2 is both (Default both).
double epd_mult(TClonesArray *epd_hits, int side=2);
bool pileup_cut(int charged_particles, int tofmatch, int tofmult, int refmult3);

int main(int argc, char **argv) {
    ROOT::EnableImplicitMT();
    // Set up QA manager, controls processes, cuts, io, etc...
    QA_Manager *manager = new QA_Manager(argc, argv);    

    
    // Trees
    TFile *outfile = new TFile(Form("%sout.root", manager->job_id.c_str()), "RECREATE");
    TTree *jet_data = new TTree("jet_data", "Jet Data");
    jet_data->SetDirectory(outfile);
    
    Event_Tree *event_tree = new Event_Tree(jet_data, "event");
    event_tree->writeable_tree();
    Jet_Tree *hardcore_jet_tree = new Jet_Tree(jet_data, "hc");
    hardcore_jet_tree->writeable_tree();
    Jet_Tree *jet_tree = new Jet_Tree(jet_data, "all");
    jet_tree->writeable_tree();
      
    
    // Set up jet finding
    Jet_Helper *jet_helper = new Jet_Helper(0.3);
    Jet_Helper *hardcore_jet_helper = new Jet_Helper(0.3);

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
        bool is_minbias = false;
        bool is_bht1_vpd30 = false;
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

        if (!manager->only_ep_finding && !is_bht1_vpd30) {   // Only run analysis on bht1_vpd30 trigger - feels weird
            continue;
        }

        // Find event plane
        if (manager->only_ep_finding && !is_minbias) {   // Only use minbias events to generate the EPD corrections - this confuses me
            continue;
        }

        if (manager->reader->centrality16() == -1) { // skip runs we can't determine centrality for
            continue;
        }

        // Pileup cut
        if (!pileup_cut(event->numberOfPrimaryTracks(), event->nBTOFMatch(), event->btofTrayMultiplicity(), event->refMult3())) {
            continue;
        }
        
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
        


        std::vector<fastjet::PseudoJet> tracks = manager->reader->pseudojets();
        std::vector<fastjet::PseudoJet> hardcore_tracks;
        for (std::vector<fastjet::PseudoJet>::iterator track = tracks.begin(); track != tracks.end(); track++) {    // Let's just play with tpc tracks for now
            jetreader::VectorInfo track_info = track->user_info<jetreader::VectorInfo>();
            if (track_info.isBemcTower()) {
                std::cout << "found tower " << track_info.towerId() << std::endl;
                
            }
            if ((!track_info.isPrimary())) {
                tracks.erase(track);
                track--;
            }
            else {
                double track_pt = track->pt();
                qa_hist.track_momentum->Fill(track_pt);
                qa_hist.track_loc->Fill(track->eta(), track->phi());
                if (track_pt > 2) {
                    hardcore_tracks.push_back(*track);
                }
            }
        }
        // std::cout << "Regular track count: " << tracks.size() << std::endl;
        // std::cout << "Hardcore track count: " << hardcore_tracks.size() << std::endl;


        std::vector<fastjet::PseudoJet> all_jets = jet_helper->find_jets(tracks);
        jet_helper->set_background_particles(tracks);
        jet_helper->fill_jet_tree(all_jets, jet_tree);
        
        std::vector<fastjet::PseudoJet> hardcore_jets = hardcore_jet_helper->find_jets(hardcore_tracks);
        hardcore_jet_helper->set_background_particles(hardcore_tracks);
        hardcore_jet_helper->fill_jet_tree(hardcore_jets, hardcore_jet_tree);
        
        
        // std::cout << "Regular jet cout: " << jet_tree->num_jets << std::endl;
        // std::cout << "Hardcore jet cout: " << hardcore_jet_tree->num_jets << std::endl;
        if (jet_tree->num_jets) {
            event_tree->fill_tree();
        }
    }

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
