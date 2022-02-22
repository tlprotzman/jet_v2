#include "setup.h"

#include "jetreader/reader/reader.h"
#include "jetreader/reader/centrality_def.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TFile.h"

#include <string>

void setup_cuts(jetreader::Reader *reader) {
    // Properties
    // Vertex Selection
    float vz_min = -50;
    float vz_max = 50;
    float vr_max = 0.5;
    std::string bad_run_list = "";

    // Track cuts
    float track_pt_max = 30;
    float track_dca_max = 3;
    float track_nhits_min = 15;
    float track_nhits_frac_min = 0.52;
    bool reject_track_on_fail = true;

    // Tower cuts
    float tower_et_max = 30;
    float tower_hadronic_correction_factor = 1;
    bool reject_tower_on_fail = true;
    bool tower_use_hadronic_correction = true;
    bool tower_use_approximate_track_tower_matching = true;
    std::string bad_tower_list = "";


    // Apply reader properties
    if (bad_run_list != "") {
        reader->eventSelector()->addBadRuns(bad_run_list);
    }

    // Vertex
    reader->eventSelector()->setVzRange(vz_min, vz_max);
    reader->eventSelector()->setVrMax(vr_max);

    // Tracks
    reader->trackSelector()->setPtMax(track_pt_max);
    reader->trackSelector()->setDcaMax(track_dca_max);
    reader->trackSelector()->setNHitsMin(track_nhits_min);
    reader->trackSelector()->setNHitsFracMin(track_nhits_frac_min);
    reader->trackSelector()->rejectEventOnPtFailure(reject_track_on_fail);

    // Tower
    if (bad_tower_list != "") {
        reader->towerSelector()->addBadTowers(bad_tower_list);
    }
    reader->useHadronicCorrection(tower_use_hadronic_correction, tower_hadronic_correction_factor);
    reader->useApproximateTrackTowerMatching(tower_use_approximate_track_tower_matching);
    reader->towerSelector()->setEtMax(tower_et_max);
    reader->towerSelector()->rejectEventOnEtFailure(reject_tower_on_fail);

    // Centrality
    reader->centrality().loadCentralityDef(jetreader::CentDefId::Run18Zr);
    reader->centrality().loadCentralityDef(jetreader::CentDefId::Run18Ru);

    reader->Init();
}

void setup_tree(TTree *tree, jet_tree_data *datum) {
    // Jet Spatial Components
    tree->Branch("jet_eta", &datum->jet_eta);
    tree->Branch("jet_phi", &datum->jet_phi);
    tree->Branch("vx", &datum->vx);
    tree->Branch("vy", &datum->vy);
    tree->Branch("vz", &datum->vz);
    tree->Branch("vpd_vz", &datum->vpd_vz);

    // Jet Momentum Components
    tree->Branch("jet_momentum", &datum->jet_momentum);
    tree->Branch("jet_momentum_median_subtracted", &datum->jet_momentum_medium_subtracted);
    tree->Branch("jet_hardcore_momentum", &datum->jet_hardcore_momentum);
    tree->Branch("hardcore_jets", &datum->hardcore_jets);
    tree->Branch("all_jets", &datum->all_jets);

    // Jet Shape Components
    tree->Branch("jet_z", &datum->jet_z);
    
    // Collision information
    tree->Branch("event_plane_east", &datum->event_plane_east);
    tree->Branch("event_plane_west", &datum->event_plane_west);
    tree->Branch("event_plane_full", &datum->event_plane_full);
    tree->Branch("centrality16", &datum->centrality);

    // System
    tree->Branch("tofmult", &datum->tofmult);
    tree->Branch("refmult3", &datum->refmult3);
    tree->Branch("trigger_id", &datum->trigger_id);
    tree->Branch("run_number", &datum->run_number);
    tree->Branch("bbc_east_rate", &datum->bbc_east_rate);
    tree->Branch("bbc_west_rate", &datum->bbc_west_rate);
}


void setup_histograms(qa_histograms *qa_hist, ep_histograms *ep_hist) {
    // Vertex info
    qa_hist->vz = new TH1D("vz", "vz", 103, -51, 51);
    qa_hist->vr = new TH2D("vr", "vr", 103,  -0.51, 0.51, 103, -0.51, 0.51);

    // System
    qa_hist->refmult3 = new TH1I("refmult3", "Ref Mult 3", 100, 0, 800);
    qa_hist->tofmult = new TH1I("tofmult", "Tof Mult", 100, 0, 1500);
    qa_hist->bbc_rate = new TH2F("bbc_rate", "BBC Rate", 100, 70000, 80000, 100, 70000, 80000);
    qa_hist->centrailty = new TH1I("centrality", "centrality", 18, -1, 16);

    // Track info
    qa_hist->track_momentum = new TH1D("track_momentum", "Track Momentum", 100, 0, 32);
    qa_hist->track_eta = new TH1D("track_eta", "Track Eta", 50, -1.1, 1.1);
    qa_hist->track_phi = new TH1D("track_phi", "Track Phi", 50, 0, TMath::TwoPi());

    // Event plane
    ep_hist->east_uncorrected = new TH1D("east_uncorrected", "east_uncorrected", 30, 0, TMath::Pi());
    ep_hist->west_uncorrected = new TH1D("west_uncorrected", "west_uncorrected", 30, 0, TMath::Pi());
    ep_hist->east_phi_corrected = new TH1D("east_phi_corrected", "east_phi_corrected", 30, 0, TMath::Pi());
    ep_hist->west_phi_corrected = new TH1D("west_phi_corrected", "west_phi_corrected", 30, 0, TMath::Pi());
    ep_hist->east_phi_psi_corrected = new TH1D("east_phi_psi_corrected", "east_phi_psi_corrected", 30, 0, TMath::Pi());
    ep_hist->west_phi_psi_corrected = new TH1D("west_phi_psi_corrected", "west_phi_psi_corrected", 30, 0, TMath::Pi());
    ep_hist->ep_correlation = new TH2D("ep_correlation", "EP Correlation", 30, 0, TMath::Pi(), 30, 0, TMath::Pi());
}

void save_histograms(qa_histograms *qa_hist, ep_histograms *ep_hist, TFile *outfile) {
    // Vertex info
    qa_hist->vz->SetDirectory(outfile);
    qa_hist->vz->Write();
    qa_hist->vr->SetDirectory(outfile);
    qa_hist->vr->Write();

    // System
    qa_hist->refmult3->SetDirectory(outfile);
    qa_hist->refmult3->Write();
    qa_hist->tofmult->SetDirectory(outfile);
    qa_hist->tofmult->Write();
    qa_hist->bbc_rate->SetDirectory(outfile);
    qa_hist->bbc_rate->Write();
    qa_hist->centrailty->SetDirectory(outfile);
    qa_hist->centrailty->Write();


    // Track info
    qa_hist->track_momentum->SetDirectory(outfile);
    qa_hist->track_momentum->Write();
    qa_hist->track_eta->SetDirectory(outfile);
    qa_hist->track_eta->Write();
    qa_hist->track_phi->SetDirectory(outfile);
    qa_hist->track_phi->Write();

    // Event plane
    ep_hist->east_uncorrected->SetDirectory(outfile);
    ep_hist->east_uncorrected->Write();
    ep_hist->west_uncorrected->SetDirectory(outfile);
    ep_hist->west_uncorrected->Write();
    ep_hist->east_phi_corrected->SetDirectory(outfile);
    ep_hist->east_phi_corrected->Write();
    ep_hist->west_phi_corrected->SetDirectory(outfile);
    ep_hist->west_phi_corrected->Write();
    ep_hist->east_phi_psi_corrected->SetDirectory(outfile);
    ep_hist->east_phi_psi_corrected->Write();
    ep_hist->west_phi_psi_corrected->SetDirectory(outfile);
    ep_hist->west_phi_psi_corrected->Write();
    ep_hist->ep_correlation->SetDirectory(outfile);
    ep_hist->ep_correlation->Write();



}