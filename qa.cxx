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

double JET_PT_CUT = 5;
double NMIP_MIN = 0.3;
double NMIP_MAX = 3;

// Calculate the sum of nMips in the EPD.  Side 0 is east, 1 is west, 2 is both (Default both).
double epd_mult(TClonesArray *epd_hits, int side=2);
bool pileup_cut(int charged_particles, int tofmatch, int tofmult, int refmult3);

int main(int argc, char **argv) {
    // Input argument parsing
    int optstring;
    int n = kMaxInt;
    bool only_ep_finding = false;
    bool nocuts = false;
    bool has_pico_list = false;
    std::string jobid = "";
    std::string picos_to_read;
    while ((optstring = getopt(argc, argv, "f:en:j:sch")) != -1) {
        switch(optstring) {
            case 'f':
                picos_to_read = std::string(optarg);
                has_pico_list = true;
                break;
            case 'e':
                only_ep_finding = true;
                break;
            case 'n':
                n = std::stoi(optarg);
                break;
            case 'j':
                jobid = std::string(optarg);
                break;
            case 's':
                picos_to_read = "myfilelist.list";
                has_pico_list = true;
                break;
            case 'c':
                nocuts = true;
                break;
            case 'h':
                std::cout << "Runs QA for Tristan's jet v2 analysis, targeting the isobar dataset\n";
                std::cout << "\t-f\tList of picos to read\n";
                std::cout << "\t-e\tOnly finds event plane, no jet finding\n";
                std::cout << "\t-n\tNumber of events to run over\n";
                std::cout << "\t-j\tJob ID to be appended to file name\n";
                std::cout << "\t-c\tRun without QA cuts\n";
                return 0;
            default: break;
        }
    }
    if (!has_pico_list) {
        std::cout << "invoke with " << argv[0] << " -f {pico list}" << std::endl;
        return -1;
    }

    std::cout << Form("Running %d events from %s%s\n\n", n, picos_to_read.c_str(), only_ep_finding ? ", just finding event plane" : "") << std::endl;
    
    // ROOT::EnableImplicitMT();

    // PARAMETERS
    int DEBUG_LEVEL = 0;
    if (getenv("JR_DEBUG_LEVEL") != nullptr) {
        DEBUG_LEVEL = atoi(getenv("JR_DEBUG_LEVEL"));
    }
    std::cout << "Running at debug level " << DEBUG_LEVEL << std::endl;
    
    std::string bad_run_list = "";


    // INITIALIZATION
    jetreader::Reader *reader = new jetreader::Reader(picos_to_read);
    setup_cuts(reader, nocuts);

    
    // Trees
    TFile *outfile = new TFile(Form("%sout.root", jobid.c_str()), "RECREATE");
    TTree *jet_data = new TTree("jet_data", "Jet Data");
    jet_data->SetDirectory(outfile);
    jet_tree_data datum;
    setup_tree(jet_data, &datum);    
    
    // Set up jet finding
    double jet_r = 0.3;
    double ghost_maxrap = 1;
    fastjet::GhostedAreaSpec area_spec(ghost_maxrap);
    fastjet::AreaDefinition jet_area(fastjet::active_area, area_spec);
    fastjet::AreaDefinition jet_area_background(fastjet::active_area_explicit_ghosts, fastjet::GhostedAreaSpec(ghost_maxrap));

    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, jet_r);
    fastjet::JetDefinition jet_subtraction_jet_def(fastjet::kt_algorithm, jet_r);

    fastjet::Selector jet_selector = fastjet::SelectorAbsRapMax(1) * (!fastjet::SelectorNHardest(2));
    fastjet::JetMedianBackgroundEstimator jet_background_estimator = fastjet::JetMedianBackgroundEstimator(jet_selector, jet_subtraction_jet_def, jet_area_background);

    fastjet::Subtractor jet_background_subtractor(&jet_background_estimator);

    jet_background_subtractor.set_use_rho_m(true);
    jet_background_subtractor.set_safe_mass(true);

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
    while (reader->next()) {
        processed_events++;
        if (processed_events % 1000 == 0) {
            std::cout << "Processed " << processed_events << " events" << std::endl;
        }
        if (processed_events > n){
            break;
        }

        StPicoEvent *event = reader->picoDst()->event();

        if (reader->centrality16() == -1) { // skip runs we can't determine centrality for
            continue;
        }

        // Event Info
        // datum.trigger_id = event->triggerIds();
        bool is_minbias = false;
        bool is_bht1_vpd30 = false;
        datum.num_triggers = 0;
        for (auto t : event->triggerIds()) {
            if (trigger_names.triggers[t] == trigger_names.minbias) {
                is_minbias = true;
            } else if (trigger_names.triggers[t] == trigger_names.bht1_vpd30) {
                is_bht1_vpd30 = true;
            }
            datum.triggers[datum.num_triggers] = t;
            datum.num_triggers++;
            if (datum.num_entries <= datum.num_triggers) {
                std::cerr << "WARNING: NOT ENOUGH SPACE FOR ALL TRIGGERS" << std::endl;
            }
        }

        if (!only_ep_finding && !is_bht1_vpd30) {   // Only run analysis on bht1_vpd30 trigger - feels weird
            continue;
        }

        // Find event plane
        if (only_ep_finding && !is_minbias) {   // Only use minbias events to generate the EPD corrections - this confuses me
            continue;
        }
        // std::cout << "Minbias: " << is_minbias << "\tbht1_vpd30: " << is_bht1_vpd30 << "\n";

        // Pileup cut
        if (!pileup_cut(event->numberOfPrimaryTracks(), event->nBTOFMatch(), event->btofTrayMultiplicity(), event->refMult3())) {
            continue;
        }
        
        datum.run_number = event->runId();

        TVector3 primary_vertex = event->primaryVertex();
        datum.vx = event->primaryVertex().X();
        datum.vy = event->primaryVertex().Y();
        datum.vz = event->primaryVertex().Z();
        datum.vpd_vz = event->vzVpd();
        datum.refmult3 = event->refMult3();
        datum.tofmatch = event->nBTOFMatch();
        datum.tofmult = event->btofTrayMultiplicity();
        datum.bbc_east_rate = event->bbcEastRate();
        datum.bbc_west_rate = event->bbcWestRate();
        datum.centrality = reader->centrality16();

        qa_hist.vz->Fill(datum.vz, event->vzVpd());
        qa_hist.vr->Fill(datum.vx, datum.vy);
        qa_hist.pileup->Fill(datum.refmult3, datum.tofmult);
        qa_hist.tofmult->Fill(datum.tofmult);
        qa_hist.bbc_rate->Fill(datum.bbc_east_rate, datum.bbc_west_rate);
        qa_hist.centrality->Fill(datum.centrality);
        
       
        TClonesArray *epd_hits = reader->picoDst()->picoArray(8);
        double nMips_sum = epd_mult(epd_hits);
        qa_hist.nMips->Fill(datum.refmult3, nMips_sum); // save nmips for jet events too?
        StEpdEpInfo ep_info = ep_finder->Results(epd_hits, primary_vertex, datum.centrality);
        ep_hist.east_uncorrected->Fill(ep_info.EastRawPsi(2));
        ep_hist.west_uncorrected->Fill(ep_info.WestRawPsi(2));
        ep_hist.east_phi_corrected->Fill(ep_info.EastPhiWeightedPsi(2));
        ep_hist.west_phi_corrected->Fill(ep_info.WestPhiWeightedPsi(2));
        ep_hist.east_phi_psi_corrected->Fill(ep_info.EastPhiWeightedAndShiftedPsi(2));
        ep_hist.west_phi_psi_corrected->Fill(ep_info.WestPhiWeightedAndShiftedPsi(2));
        ep_hist.ep_correlation->Fill(ep_info.EastPhiWeightedAndShiftedPsi(2), ep_info.WestPhiWeightedAndShiftedPsi(2));
        ep_hist.epd_resolution->Fill(ep_info.EastPhiWeightedAndShiftedPsi(2) - ep_info.WestPhiWeightedAndShiftedPsi(2), reader->centrality16());
        datum.event_plane_east = ep_info.EastPhiWeightedAndShiftedPsi(2);
        datum.event_plane_west = ep_info.WestPhiWeightedAndShiftedPsi(2);
        datum.event_plane_full = ep_info.FullPhiWeightedAndShiftedPsi(2);
        if (only_ep_finding) {
            jet_data->Fill();
            continue;
        }
        


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
                qa_hist.track_momentum->Fill(track_pt);
                qa_hist.track_loc->Fill(track->eta(), track->phi());
                if (track_pt > 2) {
                    hardcore_tracks.push_back(*track);
                }
            }
        }

        fastjet::ClusterSequenceArea hc_cs = fastjet::ClusterSequenceArea(hardcore_tracks, jet_def, jet_area);
        fastjet::ClusterSequenceArea cs = fastjet::ClusterSequenceArea(tracks, jet_def, jet_area);

        std::vector<fastjet::PseudoJet> hardcore_jets = fastjet::sorted_by_pt(hc_cs.inclusive_jets());
        std::vector<fastjet::PseudoJet> all_jets = fastjet::sorted_by_pt(cs.inclusive_jets());
        
        jet_background_estimator.set_particles(tracks);
        
        clear_vectors(&datum);
        for (fastjet::PseudoJet jet : hardcore_jets) {  
            if (jet.pt() < JET_PT_CUT) {
                break;
            }
            datum.hardcore_jets_phi[datum.num_hardcore_jets] = jet.phi();
            datum.hardcore_jets_eta[datum.num_hardcore_jets] = jet.eta();
            datum.hardcore_jets_pt[datum.num_hardcore_jets] = jet.pt();
            datum.hardcore_jets_subtracted_pt[datum.num_all_jets] = jet.pt() - jet_background_estimator.rho() * jet.area_4vector().pt();
            datum.hardcore_jets_E[datum.num_hardcore_jets] = jet.E();
            datum.num_hardcore_jets++;
            if (datum.num_entries <= datum.num_hardcore_jets) {
                break;
            }
        }

        for (fastjet::PseudoJet jet : all_jets) {
            qa_hist.jet_loc->Fill(jet.eta(), jet.phi());
            qa_hist.jet_pt_spectra->Fill(jet.pt());
            double subtracted_pt = jet.pt() - jet_background_estimator.rho() * jet.area_4vector().pt();
            qa_hist.jet_subtracted_pt_spectra->Fill(subtracted_pt);
            if (subtracted_pt < JET_PT_CUT) {
                continue;
            }
            if (jet.constituents().size() < 2) {
                continue;
            }
            if (datum.num_entries <= datum.num_all_jets) {
                continue;
            }
            datum.all_jets_phi[datum.num_all_jets] = jet.phi();
            datum.all_jets_eta[datum.num_all_jets] = jet.eta();
            datum.all_jets_pt[datum.num_all_jets] = jet.pt();
            datum.all_jets_subtracted_pt[datum.num_all_jets] = jet.pt() - jet_background_estimator.rho() * jet.area_4vector().pt();
            datum.all_jets_E[datum.num_all_jets] = jet.E();
            datum.all_jets_constituents[datum.num_all_jets] = jet.constituents().size();
            double max_pt = 0;
            for (auto constituent : jet.constituents()) {
                max_pt = constituent.pt() > max_pt ? constituent.pt() : max_pt;
            }
            datum.all_jets_z[datum.num_all_jets] = max_pt / jet.pt();
            datum.num_all_jets++;
        }
        jet_data->Fill();
        
    }
    std::cout << "Count: " << processed_events << std::endl;
    ep_finder->Finish();


    // Write data to file
    outfile->cd();
    jet_data->Write();
    save_histograms(&qa_hist, &ep_hist, outfile);
    outfile->Close();
    // delete jet_data;
    cleanup(&datum);
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
