#include "setup.h"

#include "jetreader/reader/reader.h"
#include "jetreader/reader/centrality_def.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TFile.h"

#include <string>

int NUM_ENTRIES = 10;

void setup_cuts(jetreader::Reader *reader, bool nocuts) {
    // Properties
    // Vertex Selection
    float vz_min = -35;
    float vz_max = 25;
    float vr_max = 2;
    std::string bad_run_list = "badrunlist.txt";

    // Track cuts
    float track_pt_max = 30;
    float track_dca_max = 1;
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


    if (!nocuts) {
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

    }
    
    // Centrality
    reader->centrality().loadCentralityDef(jetreader::CentDefId::Run18Zr);
    reader->centrality().loadCentralityDef(jetreader::CentDefId::Run18Ru);
    
    reader->Init();
}

void setup_tree(TTree *tree, jet_tree_data *datum) {
    // Initialize vectors
    datum->num_entries = NUM_ENTRIES;
    datum->num_hardcore_jets = 0; // Always initialize your variables... 
    datum->hardcore_jets_pt = (double*) malloc(NUM_ENTRIES * sizeof(double));
    datum->hardcore_jets_eta = (double*) malloc(NUM_ENTRIES * sizeof(double));
    datum->hardcore_jets_phi = (double*) malloc(NUM_ENTRIES * sizeof(double));
    datum->hardcore_jets_E = (double*) malloc(NUM_ENTRIES * sizeof(double));
    datum->hardcore_jets_subtracted_pt = (double*) malloc(NUM_ENTRIES * sizeof(double));
    datum->num_all_jets = 0;  // grr...
    datum->all_jets_pt = (double*) malloc(NUM_ENTRIES * sizeof(double));
    datum->all_jets_eta = (double*) malloc(NUM_ENTRIES * sizeof(double));
    datum->all_jets_phi = (double*) malloc(NUM_ENTRIES * sizeof(double));
    datum->all_jets_E = (double*) malloc(NUM_ENTRIES * sizeof(double));
    datum->all_jets_subtracted_pt = (double*) malloc(NUM_ENTRIES * sizeof(double));
    datum->all_jets_z = (double*) malloc(NUM_ENTRIES * sizeof(double));
    datum->all_jets_constituents = (UInt_t*) malloc(NUM_ENTRIES * sizeof(UInt_t));
    datum->triggers = (long*) malloc(NUM_ENTRIES * sizeof(long));

    // Vertex Components
    tree->Branch("vx", &datum->vx);
    tree->Branch("vy", &datum->vy);
    tree->Branch("vz", &datum->vz);
    tree->Branch("vpd_vz", &datum->vpd_vz);

    // Jet Components
    tree->Branch("hardcore_jets_num", &datum->num_hardcore_jets);
    tree->Branch("hardcore_jets_pt", datum->hardcore_jets_pt, "hardcore_jets_pt[hardcore_jets_num]/D");
    tree->Branch("hardcore_jets_eta", datum->hardcore_jets_eta, "hardcore_jets_eta[hardcore_jets_num]/D");
    tree->Branch("hardcore_jets_phi", datum->hardcore_jets_phi, "hardcore_jets_phi[hardcore_jets_num]/D");
    tree->Branch("hardcore_jets_E", datum->hardcore_jets_E, "hardcore_jets_E[hardcore_jets_num]/D");
    tree->Branch("hardcore_jets_subtracted_pt", datum->hardcore_jets_subtracted_pt, "hardcore_jets_subtracted_pt[hardcore_jets_num]/D");
    tree->Branch("all_jets_num", &datum->num_all_jets);
    tree->Branch("all_jets_pt", datum->all_jets_pt, "all_jets_pt[all_jets_num]/D");
    tree->Branch("all_jets_eta", datum->all_jets_eta, "all_jets_eta[all_jets_num]/D");
    tree->Branch("all_jets_phi", datum->all_jets_phi, "all_jets_phi[all_jets_num]/D");
    tree->Branch("all_jets_E", datum->all_jets_E, "all_jets_E[all_jets_num]/D");
    tree->Branch("all_jets_subtracted_pt", datum->all_jets_subtracted_pt, "all_jets_subtracted_pt[all_jets_num]/D");
    tree->Branch("all_jets_z", datum->all_jets_z, "all_jets_z[all_jets_num]/D");
    tree->Branch("all_jets_constituents", datum->all_jets_constituents, "all_jets_constituents[all_jets_num]/i");

    // Collision information
    tree->Branch("event_plane_east", &datum->event_plane_east);
    tree->Branch("event_plane_west", &datum->event_plane_west);
    tree->Branch("event_plane_full", &datum->event_plane_full);
    tree->Branch("centrality16", &datum->centrality);

    // System
    tree->Branch("tofmult", &datum->tofmult);
    tree->Branch("refmult3", &datum->refmult3);
    tree->Branch("num_triggers", &datum->num_triggers);
    tree->Branch("triggers", datum->triggers, "triggers[num_triggers]/L");
    tree->Branch("run_number", &datum->run_number);
    tree->Branch("bbc_east_rate", &datum->bbc_east_rate);
    tree->Branch("bbc_west_rate", &datum->bbc_west_rate);
}

// There must be a good way to semi-automatically generate this list...
void read_tree(TTree *tree, jet_tree_data *datum) {
    // Initialize vectors
    datum->num_entries = NUM_ENTRIES;
    datum->hardcore_jets_pt = (double*) malloc(NUM_ENTRIES * sizeof(double));
    datum->hardcore_jets_eta = (double*) malloc(NUM_ENTRIES * sizeof(double));
    datum->hardcore_jets_phi = (double*) malloc(NUM_ENTRIES * sizeof(double));
    datum->hardcore_jets_E = (double*) malloc(NUM_ENTRIES * sizeof(double));
    datum->hardcore_jets_subtracted_pt = (double*) malloc(NUM_ENTRIES * sizeof(double));
    datum->all_jets_pt = (double*) malloc(NUM_ENTRIES * sizeof(double));
    datum->all_jets_eta = (double*) malloc(NUM_ENTRIES * sizeof(double));
    datum->all_jets_phi = (double*) malloc(NUM_ENTRIES * sizeof(double));
    datum->all_jets_E = (double*) malloc(NUM_ENTRIES * sizeof(double));
    datum->all_jets_subtracted_pt = (double*) malloc(NUM_ENTRIES * sizeof(double));
    datum->all_jets_z = (double*) malloc(NUM_ENTRIES * sizeof(double));
    datum->all_jets_constituents = (UInt_t*) malloc(NUM_ENTRIES * sizeof(UInt_t));

    // Vertex Components
    tree->SetBranchAddress("vx", &datum->vx);
    tree->SetBranchAddress("vy", &datum->vy);
    tree->SetBranchAddress("vz", &datum->vz);
    tree->SetBranchAddress("vpd_vz", &datum->vpd_vz);

    // Jet Components
    tree->SetBranchAddress("hardcore_jets_num", &datum->num_hardcore_jets);
    tree->SetBranchAddress("hardcore_jets_pt", datum->hardcore_jets_pt);
    tree->SetBranchAddress("hardcore_jets_eta", datum->hardcore_jets_eta);
    tree->SetBranchAddress("hardcore_jets_phi", datum->hardcore_jets_phi);
    tree->SetBranchAddress("hardcore_jets_E", datum->hardcore_jets_E);
    tree->SetBranchAddress("hardcore_jets_subtracted_pt", datum->hardcore_jets_subtracted_pt);
    tree->SetBranchAddress("all_jets_num", &datum->num_all_jets);
    tree->SetBranchAddress("all_jets_pt", datum->all_jets_pt);
    tree->SetBranchAddress("all_jets_eta", datum->all_jets_eta);
    tree->SetBranchAddress("all_jets_phi", datum->all_jets_phi);
    tree->SetBranchAddress("all_jets_E", datum->all_jets_E);
    tree->SetBranchAddress("all_jets_subtracted_pt", datum->all_jets_subtracted_pt);
    tree->SetBranchAddress("all_jets_z", datum->all_jets_z);
    tree->SetBranchAddress("all_jets_constituents", datum->all_jets_constituents);

    // Collision information
    tree->SetBranchAddress("event_plane_east", &datum->event_plane_east);
    tree->SetBranchAddress("event_plane_west", &datum->event_plane_west);
    tree->SetBranchAddress("event_plane_full", &datum->event_plane_full);
    tree->SetBranchAddress("centrality16", &datum->centrality);

    std::vector<std::seed_seq::result_type> *test;

    // System
    tree->SetBranchAddress("tofmult", &datum->tofmult);
    tree->SetBranchAddress("refmult3", &datum->refmult3);
    // tree->SetBranchAddress("trigger_id", &test);
    tree->SetBranchAddress("run_number", &datum->run_number);
    tree->SetBranchAddress("bbc_east_rate", &datum->bbc_east_rate);
    tree->SetBranchAddress("bbc_west_rate", &datum->bbc_west_rate);
}

void clear_vectors(jet_tree_data *datum) {
    datum->num_hardcore_jets = 0;
    datum->num_all_jets = 0;
}


void setup_histograms(qa_histograms *qa_hist, ep_histograms *ep_hist) {
    // Vertex info
    qa_hist->vz = new TH2D("vz", "vz", 103, -51, 51, 103, -51, 51);
    qa_hist->vz->GetXaxis()->SetTitle("vz");
    qa_hist->vz->GetYaxis()->SetTitle("vpd vz");
    qa_hist->vr = new TH2D("vr", "vr", 103,  -0.51, 0.51, 103, -0.51, 0.51);
    qa_hist->vr->GetXaxis()->SetTitle("vx");
    qa_hist->vr->GetYaxis()->SetTitle("vy");

    // System
    qa_hist->pileup = new TH2I("pileup", "Pileup Rejection", 200, 0, 800, 200, 0, 800);
    qa_hist->pileup->GetXaxis()->SetTitle("Refmult3");
    qa_hist->pileup->GetYaxis()->SetTitle("Tofmatch");
    qa_hist->tofmult = new TH1I("tofmult", "Tof Mult", 100, 0, 1500);
    qa_hist->bbc_rate = new TH2F("bbc_rate", "BBC Rate", 100, 70000, 80000, 100, 70000, 80000);
    qa_hist->bbc_rate->GetXaxis()->SetTitle("BBC East");
    qa_hist->bbc_rate->GetYaxis()->SetTitle("BBC West");
    qa_hist->centrality = new TH1I("centrality", "centrality", 18, -1, 16);
    qa_hist->nMips = new TH2D("nMips", "nMips vs RefMult", 200, 0, 800, 200, 0, 2500);
    qa_hist->nMips->GetXaxis()->SetTitle("RefMult3");
    qa_hist->nMips->GetYaxis()->SetTitle("EPD nMips Sum");


    // Track info
    qa_hist->track_momentum = new TH1D("track_momentum", "Track Momentum", 100, 0, 32);
    qa_hist->track_loc = new TH2D("track_loc", "Track Loc", 100, -1.1, 1.1, 96, 0, TMath::TwoPi());

    // Jet info
    qa_hist->jet_pt_spectra = new TH1D("jet_momentum", "Jet Momentum", 90, 0, 45);
    qa_hist->jet_subtracted_pt_spectra = new TH1D("jet_subtracted_momentum", "Jet Subtracted Momentum", 100, -10, 45);
    qa_hist->jet_loc = new TH2D("jet_loc", "Jet Loc", 100, -1.1, 1.1, 96, 0, TMath::TwoPi());

    // Event plane
    ep_hist->east_uncorrected = new TH1D("east_uncorrected", "east_uncorrected", 30, 0, TMath::Pi());
    ep_hist->west_uncorrected = new TH1D("west_uncorrected", "west_uncorrected", 30, 0, TMath::Pi());
    ep_hist->east_phi_corrected = new TH1D("east_phi_corrected", "east_phi_corrected", 30, 0, TMath::Pi());
    ep_hist->west_phi_corrected = new TH1D("west_phi_corrected", "west_phi_corrected", 30, 0, TMath::Pi());
    ep_hist->east_phi_psi_corrected = new TH1D("east_phi_psi_corrected", "east_phi_psi_corrected", 30, 0, TMath::Pi());
    ep_hist->west_phi_psi_corrected = new TH1D("west_phi_psi_corrected", "west_phi_psi_corrected", 30, 0, TMath::Pi());
    ep_hist->ep_correlation = new TH2D("ep_correlation", "EP Correlation", 30, 0, TMath::Pi(), 30, 0, TMath::Pi());
    ep_hist->epd_resolution = new TH2D("ep_resolution", "ep_resolution", 100, -1 * TMath::Pi(), TMath::Pi(), 18,-1, 16);
    ep_hist->epd_resolution->GetXaxis()->SetTitle("#Psi_{ep}^{east}-#Psi_{ep}^{west}");
    ep_hist->epd_resolution->GetYaxis()->SetTitle("Centrality Bin");
}

void save_histograms(qa_histograms *qa_hist, ep_histograms *ep_hist, TFile *outfile) {
    // Vertex info
    qa_hist->vz->SetDirectory(outfile);
    qa_hist->vz->Write();
    qa_hist->vr->SetDirectory(outfile);
    qa_hist->vr->Write();

    // System
    qa_hist->pileup->SetDirectory(outfile);
    qa_hist->pileup->Write();
    qa_hist->tofmult->SetDirectory(outfile);
    qa_hist->tofmult->Write();
    qa_hist->bbc_rate->SetDirectory(outfile);
    qa_hist->bbc_rate->Write();
    qa_hist->centrality->SetDirectory(outfile);
    qa_hist->centrality->Write();
    qa_hist->nMips->SetDirectory(outfile);
    qa_hist->nMips->Write();


    // Track info
    qa_hist->track_momentum->SetDirectory(outfile);
    qa_hist->track_momentum->Write();
    qa_hist->track_loc->SetDirectory(outfile);
    qa_hist->track_loc->Write();
    qa_hist->jet_pt_spectra->SetDirectory(outfile);
    qa_hist->jet_pt_spectra->Write();
    qa_hist->jet_subtracted_pt_spectra->SetDirectory(outfile);
    qa_hist->jet_subtracted_pt_spectra->Write();
    qa_hist->jet_loc->SetDirectory(outfile);
    qa_hist->jet_loc->Write();



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
    ep_hist->epd_resolution->SetDirectory(outfile);
    ep_hist->epd_resolution->Write();
}

void cleanup(jet_tree_data *datum) {
    free(datum->hardcore_jets_pt);
    free(datum->hardcore_jets_eta);
    free(datum->hardcore_jets_phi);
    free(datum->hardcore_jets_E);
    free(datum->hardcore_jets_subtracted_pt);
    free(datum->all_jets_pt);
    free(datum->all_jets_eta);
    free(datum->all_jets_phi);
    free(datum->all_jets_E);
    free(datum->all_jets_subtracted_pt);
    free(datum->all_jets_z);
    free(datum->all_jets_constituents);
}