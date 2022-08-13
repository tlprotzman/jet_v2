#include "jet_helper.h"

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

double JET_PT_CUT = 3;

Jet_Helper::Jet_Helper(const float _jet_resolution) {
    this->jet_reso = _jet_resolution;

    this->area_spec = fastjet::GhostedAreaSpec(ghost_maxrap);
    this->background_area_spec = fastjet::GhostedAreaSpec(ghost_maxrap);
    this->jet_area = fastjet::AreaDefinition(fastjet::active_area, this->area_spec);
    this->background_jet_area = fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts, this->background_area_spec);

    this->jet_def = fastjet::JetDefinition(fastjet::antikt_algorithm, this->jet_reso);
    this->background_subtraction_def = fastjet::JetDefinition(fastjet::kt_algorithm, this->jet_reso);

    this->jet_selector = fastjet::SelectorAbsRapMax(1) * (!fastjet::SelectorNHardest(2)); // wtf?
    this->jet_background_estimator = new fastjet::JetMedianBackgroundEstimator(this->jet_selector, this->background_subtraction_def, this->background_jet_area);
    
    
    this->jet_background_subtractor = fastjet::Subtractor(jet_background_estimator);
    this->jet_background_subtractor.set_use_rho_m(true);
    this->jet_background_subtractor.set_safe_mass(true);
    this->cluster_sequence = nullptr;
}

Jet_Helper::~Jet_Helper() {
    if (this->cluster_sequence != nullptr) {
        delete this->cluster_sequence;
    }
    delete this->jet_background_estimator;
}

void Jet_Helper::set_background_particles(std::vector<fastjet::PseudoJet> &tracks) {
    this->jet_background_estimator->set_particles(tracks);
}

std::vector<fastjet::PseudoJet> Jet_Helper::find_jets(std::vector<fastjet::PseudoJet> &tracks) {
    if (this->cluster_sequence != nullptr) {
        delete this->cluster_sequence;
    }
    this->cluster_sequence = new fastjet::ClusterSequenceArea(tracks, this->jet_def, this->jet_area);
    return fastjet::sorted_by_pt(cluster_sequence->inclusive_jets());

}

void Jet_Helper::fill_jet_tree(std::vector<fastjet::PseudoJet> &all_jets, Jet_Tree *jet_tree) {
    jet_tree->clear_vectors();
    double background = this->jet_background_estimator->rho();
    double total_pt = 0;
    for (fastjet::PseudoJet jet : all_jets) {
        double subtracted_pt = jet.pt() - background * jet.area_4vector().pt();
        if (subtracted_pt < JET_PT_CUT) {
            continue;
        }
        if (jet.constituents().size() < 2) {
            continue;
        }
        if (jet_tree->num_entries < jet_tree->num_jets) {
            continue;
        }
        total_pt += jet.pt();
        jet_tree->jet_phi[jet_tree->num_jets] = jet.phi();
        jet_tree->jet_eta[jet_tree->num_jets] = jet.eta();
        jet_tree->jet_pt[jet_tree->num_jets] = jet.pt();
        jet_tree->jet_pt_median_subtracted[jet_tree->num_jets] = subtracted_pt;
        jet_tree->jet_E[jet_tree->num_jets] = jet.E();
        jet_tree->jet_area_pt[jet_tree->num_jets] = jet.area_4vector().pt();
        jet_tree->rho = background;
        jet_tree->num_constituents[jet_tree->num_jets] = jet.constituents().size();
        double max_pt = 0;
        for (auto constituent : jet.constituents()) {
            max_pt = constituent.pt() > max_pt ? constituent.pt() : max_pt;
        }
        jet_tree->jet_charged_z[jet_tree->num_jets] = max_pt / jet.pt();
        jet_tree->num_jets++;
    }
    // std::cout << "found "<< jet_tree->num_jets << " jets" << std::endl;
    // std::cout << "Total pt: " << total_pt << std::endl;
}

void Jet_Helper::fill_jet_tree_particles(std::vector<fastjet::PseudoJet> &all_jets, Jet_Tree *jet_tree) {
    jet_tree->clear_vectors();
    double total_pt = 0;
    for (fastjet::PseudoJet jet : all_jets) {
        if (jet_tree->num_entries < jet_tree->num_jets) {
            continue;
        }
        total_pt += jet.pt();
        jet_tree->jet_phi[jet_tree->num_jets] = jet.phi();
        jet_tree->jet_eta[jet_tree->num_jets] = jet.eta();
        jet_tree->jet_pt[jet_tree->num_jets] = jet.pt();
        jet_tree->jet_pt_median_subtracted[jet_tree->num_jets] = jet.pt();
        jet_tree->jet_E[jet_tree->num_jets] = 0;
        jet_tree->jet_area_pt[jet_tree->num_jets] = 0;
        jet_tree->rho = 0;
        jet_tree->num_constituents[jet_tree->num_jets] = 0;
        jet_tree->jet_charged_z[jet_tree->num_jets] = 0;
        jet_tree->num_jets++;
    }
    // std::cout << "found "<< jet_tree->num_jets << " jets" << std::endl;
    // std::cout << "Total pt: " << total_pt << std::endl;
}