#ifndef SETUP_H
#define SETUP_H

#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"

#include <vector>

#include "fastjet/PseudoJet.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"


#include "jetreader/reader/reader.h"

/*
 Centrality binning:
 Bin       Centrality (16)   Centrality (9)
16           80-100%
 0           75- 80%            70-80%
 1           70- 75%            60-70%
 2           65- 70%            50-60%
 3           60- 65%            40-50%
 4           55- 60%            30-40%
 5           50- 55%            20-30%
 6           45- 50%            10-20%
 7           40- 45%             5-10%
 8           35- 40%             0- 5%
 9           30- 35%
10           25- 30%
11           20- 25%
12           15- 20%
13           10- 15%
14            5- 10%
15            0-  5%
*/

// Tree variables

typedef struct {
    double jet_eta, jet_phi;
    double jet_hardcore_momentum, jet_momentum, jet_momentum_medium_subtracted;
    // std::vector<double> *hardcore_jets_pt;
    // std::vector<double> *hardcore_jets_eta;
    // std::vector<double> *hardcore_jets_phi;
    // std::vector<double> *hardcore_jets_E;
    // std::vector<double> *hardcore_jets_subtracted_pt;
    // std::vector<double> *hardcore_jets_z;
    // std::vector<size_t> *hardcore_jets_constituents;
    // std::vector<double> *all_jets_pt;
    // std::vector<double> *all_jets_eta;
    // std::vector<double> *all_jets_phi;
    // std::vector<double> *all_jets_E;
    // std::vector<double> *all_jets_subtracted_pt;
    // std::vector<double> *all_jets_z;
    // std::vector<size_t> *all_jets_constituents;
    UInt_t num_hardcore_jets, num_all_jets;
    double* hardcore_jets_pt;
    double* hardcore_jets_eta;
    double* hardcore_jets_phi;
    double* hardcore_jets_E;
    double* hardcore_jets_subtracted_pt;
    double* hardcore_jets_z;
    UInt_t* hardcore_jets_constituents;
    double* all_jets_pt;
    double* all_jets_eta;
    double* all_jets_phi;
    double* all_jets_E;
    double* all_jets_subtracted_pt;
    double* all_jets_z;
    UInt_t* all_jets_constituents;
    int num_entries;

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
void read_tree(TTree *tree, jet_tree_data *datum);
void clear_vectors(jet_tree_data *datum);
void setup_histograms(qa_histograms *qa_hist, ep_histograms *ep_hist);
void save_histograms(qa_histograms *qa_hist, ep_histograms *ep_hist, TFile *outfile);
void cleanup(jet_tree_data *datum);


#endif // SETUP_H
