#ifndef JET_HELPER_H
#define JET_HELPER_H
// Probably don't need all of these, don't really feel like figuring that out though...
#include <fastjet/AreaDefinition.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/tools/BackgroundEstimatorBase.hh>
#include <fastjet/tools/GridMedianBackgroundEstimator.hh>
#include <fastjet/tools/JetMedianBackgroundEstimator.hh>
#include <fastjet/tools/Subtractor.hh>

#include <vector>

#include "jet_tree.h"


class Jet_Helper {
    private: 
        const float ghost_maxrap = 1;
        fastjet::GhostedAreaSpec area_spec;
        fastjet::GhostedAreaSpec background_area_spec;
        fastjet::AreaDefinition jet_area;
        fastjet::AreaDefinition background_jet_area;

        fastjet::JetDefinition jet_def;
        fastjet::JetDefinition background_subtraction_def;

        fastjet::Selector jet_selector;
        fastjet::JetMedianBackgroundEstimator *jet_background_estimator;
        fastjet::JetMedianBackgroundEstimator hmm;

        fastjet::Subtractor jet_background_subtractor;
        fastjet::ClusterSequenceArea *cluster_sequence;

    public: 
        float jet_reso;
        
        Jet_Helper(const float _jet_resolution);
        ~Jet_Helper();
        void set_background_particles(std::vector<fastjet::PseudoJet> &tracks);
        std::vector<fastjet::PseudoJet> find_jets(std::vector<fastjet::PseudoJet> &tracks);
        void fill_jet_tree_particles(std::vector<fastjet::PseudoJet> &all_jets, Jet_Tree *jet_tree);
        void fill_jet_tree(std::vector<fastjet::PseudoJet> &all_jets, Jet_Tree *jet_tree);

};

#endif // JET_HELPER_H