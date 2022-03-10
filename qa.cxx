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

// #pragma link C++ class std::vector<fastjet::PseudoJet>+;


int main(int argc, char **argv) {
    // Build filelist
    // int num_input_files = std::stoi(getenv("INPUTFILECOUNT"));
    // std::ofstream myfilelist("myfilelist.list", std::ofstream::out);
    // for (size_t i = 0; i < num_input_files; i++) {
    //     myfilelist << getenv(Form("INPUTFILE%d", i)) << std::endl;
    //     std::cout << getenv(Form("INPUTFILE%d", i)) << std::endl;
    // }
    // myfilelist.close();
    // Input argument parsing
    if (argc < 3) {
        std::cout << "invoke with " << argv[0] << " -f {pico list}" << std::endl;
        return -1;
    }
    int optstring;
    int n = kMaxInt;
    bool only_ep_finding = false;
    std::string jobid = "";
    std::string picos_to_read;
    while ((optstring = getopt(argc, argv, "f:en:j:s")) != -1) {
        switch(optstring) {
            case 'f':
                picos_to_read = std::string(optarg);
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
                break;
            default: break;
        }
    }
    std::cout << Form("Running %d events from %s%s\n\n", n, picos_to_read.c_str(), only_ep_finding ? ", just finding event plane" : "") << std::endl;
    // std::ifstream filelist(picos_to_read);
    // std::string line;
    // while (filelist >> line) {
    //     std::cout << line << std::endl;
    // }
    // filelist.close();
    
    ROOT::EnableImplicitMT();

    // PARAMETERS
    int DEBUG_LEVEL = 0;
    if (getenv("JR_DEBUG_LEVEL") != nullptr) {
        DEBUG_LEVEL = atoi(getenv("JR_DEBUG_LEVEL"));
    }
    std::cout << "Running at debug level " << DEBUG_LEVEL << std::endl;
    
    // std::string picos_to_read = "test.list"; 
    // std::string picos_to_read = "smalltest.list"; 
    // std::string picos_to_read = "/data/star/production_isobar_2018/ReversedFullField/P20ic/2018/083/19083049/st_physics_19083049_raw_1000011.picoDst.root";
    std::string bad_run_list = "";


    // INITIALIZATION
    jetreader::Reader *reader = new jetreader::Reader(picos_to_read);
    setup_cuts(reader);
    
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
        if (processed_events > n){
            break;
        }

        StPicoEvent *event = reader->picoDst()->event();

        if (reader->centrality16() == -1) { // skip runs we can't determine centrality for
            continue;
        }

        // Event Info
        datum.trigger_id = event->triggerIds();
        datum.run_number = event->runId();

        TVector3 primary_vertex = event->primaryVertex();
        datum.vx = event->primaryVertex().X();
        datum.vy = event->primaryVertex().Y();
        datum.vz = event->primaryVertex().Z();
        qa_hist.vz->Fill(datum.vz);
        qa_hist.vr->Fill(datum.vx, datum.vy);
        datum.vpd_vz = event->vzVpd();

        datum.refmult3 = event->refMult3();
        qa_hist.refmult3->Fill(datum.refmult3);
        datum.tofmult = event->btofTrayMultiplicity();
        qa_hist.tofmult->Fill(datum.tofmult);
        datum.bbc_east_rate = event->bbcEastRate();
        datum.bbc_west_rate = event->bbcWestRate();
        qa_hist.bbc_rate->Fill(datum.bbc_east_rate, datum.bbc_west_rate);
        datum.centrality = reader->centrality16();
        qa_hist.centrailty->Fill(datum.centrality);

        
        // Find event plane
        TClonesArray *epd_hits = reader->picoDst()->picoArray(8);
        StEpdEpInfo ep_info = ep_finder->Results(epd_hits, primary_vertex, datum.centrality);
        ep_hist.east_uncorrected->Fill(ep_info.EastRawPsi(2));
        ep_hist.west_uncorrected->Fill(ep_info.WestRawPsi(2));
        ep_hist.east_phi_corrected->Fill(ep_info.EastPhiWeightedPsi(2));
        ep_hist.west_phi_corrected->Fill(ep_info.WestPhiWeightedPsi(2));
        ep_hist.east_phi_psi_corrected->Fill(ep_info.EastPhiWeightedAndShiftedPsi(2));
        ep_hist.west_phi_psi_corrected->Fill(ep_info.WestPhiWeightedAndShiftedPsi(2));
        ep_hist.ep_correlation->Fill(ep_info.EastPhiWeightedAndShiftedPsi(2), ep_info.WestPhiWeightedAndShiftedPsi(2));
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
                qa_hist.track_eta->Fill(track->eta());
                qa_hist.track_phi->Fill(track->phi());
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
            // std::cout << "hardcore: " << jet.pt() << std::endl;          
            datum.hardcore_jets_phi[datum.num_hardcore_jets] = jet.phi();
            datum.hardcore_jets_eta[datum.num_hardcore_jets] = jet.eta();
            datum.hardcore_jets_pt[datum.num_hardcore_jets] = jet.pt();
            datum.hardcore_jets_E[datum.num_hardcore_jets] = jet.E();
            // std::cout << datum.num_hardcore_jets << "\t" << datum.hardcore_jets_pt[datum.num_hardcore_jets] << std::endl;
            // datum.hardcore_jets_constituents->push_back(jet.constituents().size());
            // std::cout << jet.constituents().size() << std::endl;
            // double max_pt = 0;
            // for (auto constituent : jet.constituents()) {
            //     max_pt = constituent.pt() > max_pt ? constituent.pt() : max_pt;
            // }
            // datum.hardcore_jets_z->push_back(max_pt / jet.pt());
            datum.num_hardcore_jets++;
            if (datum.num_entries <= datum.num_hardcore_jets) {
                break;
            }
        }

        for (fastjet::PseudoJet jet : all_jets) {
            if (jet.pt() < 1) {
                break;
            }
            if (jet.constituents().size() < 2) {
                continue;
            }
            if (datum.num_entries <= datum.num_all_jets) {
                // std::cout << "overflow" << std::endl;
                break;
            }
            // std::cout << "all: " << jet.pt() << std::endl;          
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
        // std::cout << "filled? \n\n";
        
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

