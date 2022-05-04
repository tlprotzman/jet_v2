#include "tree_manager.h"

class Jet_Tree: public Tree_Manager {
    public:
        UInt_t num_jets;
        double *jet_pt;
        double *jet_pt_median_subtracted;
        double *jet_eta;
        double *jet_phi;
        double *jet_E;
        double *jet_charged_z;
        double *jet_neutral_z;
        double *jet_neutral_fraction;
        double *jet_area;
        double *rho;

    Jet_Tree(TTree *_tree, std::string _tag);
    ~Jet_Tree();

    int writable_tree();
    int readable_tree();
    void clear_vectors();  
};