#include "setup.h"

#include "jetreader/reader/reader.h"
#include "jetreader/reader/centrality_def.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TFile.h"
#include "TROOT.h"

#include <stdexcept>
#include <string>

int NUM_ENTRIES = 10;

// QA_Manager initializes the QA code to run over the data in the desired way
QA_Manager::QA_Manager(int argc, char** argv) {
    this->max_events = 0;
    this->max_events_set = false;
    this->do_cuts = true;
    this->has_pico_list = false;
    this->only_ep_finding = false;
    this->job_id = "";
    this->picos_to_read = "";
    this->debug_level = 0;

    if (this->setup_tasks(argc, argv)) {    // Set up tasks
        throw std::runtime_error("Error initializing QA_Manager");
    }
    this->reader = new jetreader::Reader(this->picos_to_read); // Set up jetreader
    this->setup_cuts();
}

QA_Manager::~QA_Manager() {
    delete this->reader;
}

// setup_tasks parses the arguments and sets appropriate variables
int QA_Manager::setup_tasks(int argc, char** argv) {
    int optstring;
    while ((optstring = getopt(argc, argv, "f:en:j:schm")) != -1) {
        switch(optstring) {
            case 'f':
                this->picos_to_read = std::string(optarg);
                this->has_pico_list = true;
                break;
            case 'e':
                this->only_ep_finding = true;
                break;
            case 'n':
                this->max_events = std::stoi(optarg);
                this->max_events_set = true;
                break;
            case 'j':
                this->job_id = std::string(optarg);
                break;
            case 's':
                picos_to_read = "myfilelist.list";
                has_pico_list = true;
                break;
            case 'c':
                this->do_cuts = false;
                break;
            case 'm':
                ROOT::EnableImplicitMT();
                break;
            case 'h':
                std::cout << "Runs QA for Tristan's jet v2 analysis, targeting the isobar dataset\n";
                std::cout << "\t-f\tList of picos to read\n";
                std::cout << "\t-e\tOnly finds event plane, no jet finding\n";
                std::cout << "\t-n\tNumber of events to run over\n";
                std::cout << "\t-j\tJob ID to be appended to file name\n";
                std::cout << "\t-c\tRun without QA cuts\n";
                return -1;
            default: break;
        }
    }
    if (!this->has_pico_list) {
        std::cout << "invoke with " << argv[0] << " -f {pico list}" << std::endl;
        return -1;
    }
    return 0;
    // std::cout << Form("Running %d events from %s%s\n\n", n, picos_to_read.c_str(), only_ep_finding ? ", just finding event plane" : "") << std::endl;
}

int QA_Manager::setup_io() {
    return 0;
}

int QA_Manager::setup_cuts() {
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


    if (this->do_cuts) {
        // Apply reader properties
        if (bad_run_list != "") {
            this->reader->eventSelector()->addBadRuns(bad_run_list);
        }

        // Vertex
        this->reader->eventSelector()->setVzRange(vz_min, vz_max);
        this->reader->eventSelector()->setVrMax(vr_max);

        // Tracks
        this->reader->trackSelector()->setPtMax(track_pt_max);
        this->reader->trackSelector()->setDcaMax(track_dca_max);
        this->reader->trackSelector()->setNHitsMin(track_nhits_min);
        this->reader->trackSelector()->setNHitsFracMin(track_nhits_frac_min);
        this->reader->trackSelector()->rejectEventOnPtFailure(reject_track_on_fail);

        // Tower
        if (bad_tower_list != "") {
            this->reader->towerSelector()->addBadTowers(bad_tower_list);
        }
        for (uint32_t i = 0; i <= 4800; i++) {   // Don't use any towers
            this->reader->towerSelector()->addBadTower(i);
        }
        this->reader->useHadronicCorrection(tower_use_hadronic_correction, tower_hadronic_correction_factor);
        this->reader->useApproximateTrackTowerMatching(tower_use_approximate_track_tower_matching);
        this->reader->towerSelector()->setEtMax(tower_et_max);
        this->reader->towerSelector()->rejectEventOnEtFailure(reject_tower_on_fail);

    }
    
    // Centrality
    this->reader->centrality().loadCentralityDef(jetreader::CentDefId::Run18Zr);
    this->reader->centrality().loadCentralityDef(jetreader::CentDefId::Run18Ru);
    
    this->reader->Init();
    return 0;
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
    qa_hist->pileup = new TH2I("pileup", "Pileup Rejection", 200, 0, 800, 200, 0, 1600);
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
    ep_hist->epd_resolution = new TH2D("ep_resolution", "ep_resolution", 100, -1 * TMath::Pi(), TMath::Pi(), 17, -0.5, 16.5);
    ep_hist->epd_resolution->GetXaxis()->SetTitle("#Psi_{ep}^{east}-#Psi_{ep}^{west}");
    ep_hist->epd_resolution->GetYaxis()->SetTitle("Centrality Bin");

    // Particle Anisotropy
    qa_hist->rho = new TH2D("eventwise_background", "Event #rho_{ch}", 250, 0, 25, 17, -0.5, 16.5);
    qa_hist->rho->GetXaxis()->SetTitle("#rho_{0}");
    qa_hist->rho->GetYaxis()->SetTitle("Centrality");
    qa_hist->v2 = new TH2D("eventwise_v2", "Event v_{2, ch}", 250, -0.5, 0.5, 17, -0.5, 16.5);
    qa_hist->v2->GetXaxis()->SetTitle("v_{2}");
    qa_hist->v2->GetYaxis()->SetTitle("Centrality");
    qa_hist->v3 = new TH2D("eventwise_v3", "Event v_{3, ch}", 250, -0.5, 0.5, 17, -0.5, 16.5);
    qa_hist->v3->GetXaxis()->SetTitle("v_{3}");
    qa_hist->v3->GetYaxis()->SetTitle("centrality");
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

    qa_hist->rho->SetDirectory(outfile);
    qa_hist->rho->Write();
    qa_hist->v2->SetDirectory(outfile);
    qa_hist->v2->Write();
    qa_hist->v3->SetDirectory(outfile);
    qa_hist->v3->Write();
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
