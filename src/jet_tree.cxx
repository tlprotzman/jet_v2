#include "jet_tree.h"

#include <TROOT.h>
#include <TTree.h>
#include <iostream>

Jet_Tree::Jet_Tree(TTree *_tree, std::string _tag) : Tree_Manager(_tree, _tag) {
    if (this->tree == nullptr) {
        throw std::runtime_error(Form("%s:%d tree is null", __FILE__, __LINE__));
    }
    this->num_constituents = (uint*)malloc(this->num_entries * sizeof(uint));
    this->jet_pt = (double*)malloc(this->num_entries * sizeof(double));
    this->jet_pt_median_subtracted = (double*)malloc(this->num_entries * sizeof(double));
    this->jet_eta = (double*)malloc(this->num_entries * sizeof(double));
    this->jet_phi = (double*)malloc(this->num_entries * sizeof(double));
    this->jet_E = (double*)malloc(this->num_entries * sizeof(double));
    this->jet_charged_z = (double*)malloc(this->num_entries * sizeof(double));
    this->jet_neutral_z = (double*)malloc(this->num_entries * sizeof(double));
    this->jet_neutral_fraction = (double*)malloc(this->num_entries * sizeof(double));
    this->jet_area_pt = (double*)malloc(this->num_entries * sizeof(double));
    this->leading_track = (double*)malloc(this->num_entries * sizeof(double));
}

Jet_Tree::~Jet_Tree() {
    free(this->num_constituents);
    free(this->jet_pt);
    free(this->jet_pt_median_subtracted);
    free(this->jet_eta);
    free(this->jet_phi);
    free(this->jet_E);
    free(this->jet_charged_z);
    free(this->jet_neutral_z);
    free(this->jet_neutral_fraction);
    free(this->jet_area_pt);
    free(this->leading_track);
}

int Jet_Tree::writeable_tree() {
    this->tree->Branch(Form("%s_num_jets", this->tag), &this->num_jets);
    this->tree->Branch(Form("%s_num_constituents", this->tag), this->num_constituents, Form("%s_num_constituents[%s_num_jets]/i", this->tag, this->tag));
    this->tree->Branch(Form("%s_jet_pt", this->tag), this->jet_pt, Form("%s_jet_pt[%s_num_jets]/D", this->tag, this->tag));
    this->tree->Branch(Form("%s_jet_pt_median_subtracted", this->tag), this->jet_pt_median_subtracted, Form("%s_jet_pt_median_subtracted[%s_num_jets]/D", this->tag, this->tag));
    this->tree->Branch(Form("%s_jet_eta", this->tag), this->jet_eta, Form("%s_jet_eta[%s_num_jets]/D", this->tag, this->tag));
    this->tree->Branch(Form("%s_jet_phi", this->tag), this->jet_phi, Form("%s_jet_phi[%s_num_jets]/D", this->tag, this->tag));
    this->tree->Branch(Form("%s_jet_E", this->tag), this->jet_E, Form("%s_jet_E[%s_num_jets]/D", this->tag, this->tag));
    this->tree->Branch(Form("%s_jet_charged_z", this->tag), this->jet_charged_z, Form("%s_jet_charged_z[%s_num_jets]/D", this->tag, this->tag));
    this->tree->Branch(Form("%s_jet_neutral_z", this->tag), this->jet_neutral_z, Form("%s_jet_neutral_z[%s_num_jets]/D", this->tag, this->tag));
    this->tree->Branch(Form("%s_jet_neutral_fraction", this->tag), this->jet_neutral_fraction, Form("%s_jet_neutral_fraction[%s_num_jets]/D", this->tag, this->tag));
    this->tree->Branch(Form("%s_jet_area_pt", this->tag), this->jet_area_pt, Form("%s_jet_area_pt[%s_num_jets]/D", this->tag, this->tag));
    this->tree->Branch(Form("%s_leading_track", this->tag), this->leading_track, Form("%s_leading_track[%s_num_jets]/D", this->tag, this->tag));
    this->tree->Branch(Form("%s_rho", this->tag), &this->rho);
    return 0;
}

int Jet_Tree::readable_tree() {
    std::cout << "loading branch with tag " << this->tag << std::endl;
    this->tree->SetBranchAddress(Form("%s_num_jets", this->tag), &this->num_jets);
    this->tree->SetBranchAddress(Form("%s_num_constituents", this->tag), this->num_constituents);
    this->tree->SetBranchAddress(Form("%s_jet_pt", this->tag), this->jet_pt);
    this->tree->SetBranchAddress(Form("%s_jet_pt_median_subtracted", this->tag), this->jet_pt_median_subtracted);
    this->tree->SetBranchAddress(Form("%s_jet_eta", this->tag), this->jet_eta);
    this->tree->SetBranchAddress(Form("%s_jet_phi", this->tag), this->jet_phi);
    std::cout << "wtf fuck you linking just workkkkk" << std::endl;
    this->tree->SetBranchAddress(Form("%s_jet_E", this->tag), this->jet_E);
    this->tree->SetBranchAddress(Form("%s_jet_charged_z", this->tag), this->jet_charged_z);
    this->tree->SetBranchAddress(Form("%s_jet_neutral_z", this->tag), this->jet_neutral_z);
    this->tree->SetBranchAddress(Form("%s_jet_neutral_fraction", this->tag), this->jet_neutral_fraction);
    this->tree->SetBranchAddress(Form("%s_jet_area_pt", this->tag), this->jet_area_pt);
    this->tree->SetBranchAddress(Form("%s_rho", this->tag), &this->rho);
    return 0;
};

void Jet_Tree::clear_vectors() {
    this->num_jets = 0;
}