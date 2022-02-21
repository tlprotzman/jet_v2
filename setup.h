#ifndef SETUP_H
#define SETUP_H

#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "jetreader/reader/reader.h"

// Tree variables
typedef struct {
    double jet_eta, jet_phi;
    double jet_hardcore_momentum, jet_momentum, jet_momentum_medium_subtracted;
    double jet_z;
    double vx, vy, vz;
    double event_plane_east, event_plane_west, event_plane_full;
    std::vector<std::seed_seq::result_type> trigger_id;
    double vpd_vz;                  // check against vz, for pileup
    unsigned short tofmult;         // for pileup cut
    int refmult3;
    int centrality;
    int run_number;
    float bbc_east_rate;
    float bbc_west_rate;
} jet_tree_data;

// QA Histograms
typedef struct {
    TH1 *vz;
    TH1 *vpd_vz;
    TH2 *vr;
    TH1 *track_momentum;
    TH1 *track_eta;
    TH1 *track_phi;
    // Calo_t calo_data; // figure out later
    TH1 *centrailty;
    TH1 *refmult3;
    TH1 *tofmult;
    TH2 *bbc_rate;
} qa_histograms;

typedef struct {
    TH1 *east_uncorrected;
    TH1 *west_uncorrected;
    TH1 *east_phi_corrected;
    TH1 *west_phi_corrected;
    TH1 *east_phi_psi_corrected;
    TH1 *west_phi_psi_corrected;
    TH2 *ep_correlation;
} ep_histograms;

void setup_cuts(jetreader::Reader *reader);
void setup_tree(TTree *tree, jet_tree_data *datum);
void setup_histograms(qa_histograms *qa_hist, ep_histograms *ep_hist);
void save_histograms(qa_histograms *qa_hist, ep_histograms *ep_hist, TFile *outfile);

#endif // SETUP_H
