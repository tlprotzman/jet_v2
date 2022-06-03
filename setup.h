#ifndef SETUP_H
#define SETUP_H

#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"

#include <vector>
#include <string>
#include <fstream>

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


class QA_Manager {
private:
    // Process selection
    bool do_cuts;
    bool has_pico_list;
    std::string picos_to_read;

    // File IO
    int debug_level;
    // std::ostream debug;

    // Cuts

    // Jet reader tools

public:
    bool max_events_set;
    int max_events;
    bool only_ep_finding;
    std::string job_id;
    TFile out_file;

    jetreader::Reader *reader;

    QA_Manager(int argc, char** argv);
    ~QA_Manager();

    int setup_tasks(int argc, char** argv);
    int setup_io();
    int setup_cuts();
    int setup_trees();
};


// Tree variables
typedef struct {
    double jet_eta, jet_phi;
    double jet_hardcore_momentum, jet_momentum, jet_momentum_medium_subtracted;
    UInt_t num_hardcore_jets, num_all_jets;
    double* hardcore_jets_pt;
    double* hardcore_jets_eta;
    double* hardcore_jets_phi;
    double* hardcore_jets_E;
    double* hardcore_jets_subtracted_pt;
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
    // std::vector<std::seed_seq::result_type> trigger_id;
    UInt_t num_triggers;
    long *triggers;
    double vpd_vz;                  // check against vz, for pileup
    unsigned short tofmult;         // for pileup cut
    unsigned short tofmatch;         // for pileup cut
    int refmult3;
    int centrality;
    int run_number;
    float bbc_east_rate;
    float bbc_west_rate;
} jet_tree_data;

// QA Histograms
typedef struct {
    TH2 *vz;
    TH2 *vr;
    TH1 *track_momentum;
    TH2 *track_loc;
    TH1 *jet_pt_spectra;
    TH1 *jet_subtracted_pt_spectra;
    TH2 *jet_loc;
    // Calo_t calo_data; // figure out later
    TH1 *centrality;
    TH2 *pileup;
    TH1 *tofmult;
    TH2 *bbc_rate;
    TH2 *nMips;
    TH2 *rho;
    TH2 *v2;
    TH2 *v3;
} qa_histograms;

typedef struct {
    TH1 *east_uncorrected;
    TH1 *west_uncorrected;
    TH1 *east_phi_corrected;
    TH1 *west_phi_corrected;
    TH1 *east_phi_psi_corrected;
    TH1 *west_phi_psi_corrected;
    TH2 *ep_correlation;
    TH2 *epd_resolution; 
} ep_histograms;

// void setup_cuts(jetreader::Reader *reader, bool nocuts);
void setup_tree(TTree *tree, jet_tree_data *datum);
void read_tree(TTree *tree, jet_tree_data *datum);
void clear_vectors(jet_tree_data *datum);
void setup_histograms(qa_histograms *qa_hist, ep_histograms *ep_hist);
void save_histograms(qa_histograms *qa_hist, ep_histograms *ep_hist, TFile *outfile);
void cleanup(jet_tree_data *datum);


#endif // SETUP_H
