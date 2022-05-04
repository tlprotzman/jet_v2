#include "jet_tree.h"

#include "TROOT.h"
#include "TTree.h"

Jet_Tree::Jet_Tree(TTree *_tree, std::string _tag) {
    Tree_Manager(_tree, _tag);
    this->jet_pt = (double*)malloc(this->num_entries * sizeof(double));
    this->jet_pt_median_subtracted = (double*)malloc(this->num_entries * sizeof(double));
    this->jet_eta = (double*)malloc(this->num_entries * sizeof(double));
    this->jet_phi = (double*)malloc(this->num_entries * sizeof(double));
    this->jet_E = (double*)malloc(this->num_entries * sizeof(double));
    this->jet_charged_z = (double*)malloc(this->num_entries * sizeof(double));
    this->jet_neutral_z = (double*)malloc(this->num_entries * sizeof(double));
    this->jet_neutral_fraction = (double*)malloc(this->num_entries * sizeof(double));
    this->jet_area = (double*)malloc(this->num_entries * sizeof(double));
    this->rho = (double*)malloc(this->num_entries * sizeof(double));
}

Jet_Tree::~Jet_Tree() {
    free(this->jet_pt);
    free(this->jet_pt_median_subtracted);
    free(this->jet_eta);
    free(this->jet_phi);
    free(this->jet_E);
    free(this->jet_charged_z);
    free(this->jet_neutral_z);
    free(this->jet_neutral_fraction);
    free(this->jet_area);
    free(this->rho);
}

int Jet_Tree::writable_tree() {
    this->tree->Branch(Form("%s/num_jets", this->tag), &this->num_jets);
    this->tree->Branch(Form("%s/jet_pt", this->tag), this->jet_pt, Form("%s/jet_pt[%s/num_jets]", this->tag, this->tag));
    this->tree->Branch(Form("%s/jet_pt_median_subtracted", this->tag), this->jet_pt_median_subtracted, Form("%s/jet_pt_median_subtracted[%s/num_jets]", this->tag, this->tag));
    this->tree->Branch(Form("%s/jet_eta", this->tag), this->jet_eta, Form("%s/jet_eta[%s/num_jets]", this->tag, this->tag));
    this->tree->Branch(Form("%s/jet_phi", this->tag), this->jet_phi, Form("%s/jet_phi[%s/num_jets]", this->tag, this->tag));
    this->tree->Branch(Form("%s/jet_E", this->tag), this->jet_E, Form("%s/jet_E[%s/num_jets]", this->tag, this->tag));
    this->tree->Branch(Form("%s/jet_charged_z", this->tag), this->jet_charged_z, Form("%s/jet_charged_z[%s/num_jets]", this->tag, this->tag));
    this->tree->Branch(Form("%s/jet_neutral_z", this->tag), this->jet_neutral_z, Form("%s/jet_neutral_z[%s/num_jets]", this->tag, this->tag));
    this->tree->Branch(Form("%s/jet_neutral_fraction", this->tag), this->jet_neutral_fraction, Form("%s/jet_neutral_fraction[%s/num_jets]", this->tag, this->tag));
    this->tree->Branch(Form("%s/jet_area", this->tag), this->jet_area, Form("%s/jet_area[%s/num_jets]", this->tag, this->tag));
    this->tree->Branch(Form("%s/rho", this->tag), this->rho, Form("%s/rho[%s/num_jets]", this->tag, this->tag));
    return 0;
}

void Jet_Tree::clear_vectors() {
    this->num_jets = 0;
}