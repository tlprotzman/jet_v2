#ifndef V2_ANALYZER
#define V2_ANALYZER

#include "event_tree.h"
#include "jet_tree.h"
#include "tree_manager.h"

#include <TROOT.h>
#include <TChain.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TGraphErrors.h>
#include <TRandom.h>

#include <string>
#include <vector>
#include <map>
#include <set>

class v2_analyzer {
public:
    v2_analyzer(double jet_resolution, std::string file_list);
    void run_analysis(int n=0);
    void save_histograms(std::string file_name);

private:
    void load_files();
    void register_histograms();
    int hardcore_matched(int index);
    void event_loop(int n=0);
    void calculate_v2();

    std::string file_list;
    double jet_reso;
    bool do_hardcore_matched;
    float hardcore_threshold;
    int centrality_lower_bound = 4;
    int centrality_upper_bound = 11;
    double hc_pt_min = 10; // GeV
    double assumed_v2 = 0.04; // Assumed charged particle v2 at mid rapidity
    
    uint num_events;
    double event_plane_resolution;

    // Definte pt binning for v2
    int n_pt_bins_v2 = 10;
    double pt_bins[11] = {0, 5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 30, 40};
    double avg_cos[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    double n_entries[8] = {0, 0, 0, 0, 0, 0, 0, 0};

    TChain *chain;
    std::vector<TH1*> *histograms;
    Event_Tree *event_tree;
    Jet_Tree *jets;
    Jet_Tree *hc_jets;
    TRandom *rng;
    std::set<int> *used_hc_jets;
    

    // Histograms
    TH2 *jet_relative_angle_no_background_mod;
    TH2 *jet_relative_angle_extra_background_mod;
    TH2 *jet_relative_angle_hc;
    TH2 *jet_relative_angle_all;
    TH2 *jet_relative_angle_hc_only;
    TH1 *jet_pt_spectra_not_matched;
    TH1 *jet_pt_spectra;
    TH1 *hc_pt_spectra;
    TH2 *ue_subtraction;
    TH3 *matching_performance;
    TH1 *rho_distribution;
    TH2 *subtraction_relative_ep;
    TH2 *ep_correlation;
    TH1 *ep_diff;
    TH2 *ep_resolution_cent_binned;

    TGraphErrors *jet_v2_no_modulation;
    TGraphErrors *jet_v2_extra_modulation;
    TGraphErrors *jet_v2_hc_only;
    TGraphErrors *jet_v2_all;
    TGraphErrors *jet_v2_hc;

};


#endif // V2_ANALYZER