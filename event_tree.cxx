#include "event_tree.h"

#include "TROOT.h"
#include "TTree.h"

#include <iostream>

Event_Tree::Event_Tree(TTree *_tree, std::string _tag) : Tree_Manager(_tree, _tag) {
    Tree_Manager(_tree, _tag);
    this->triggers = (long*)malloc(this->num_entries * sizeof(long));
}

Event_Tree::~Event_Tree() {
    free(this->triggers);
}

int Event_Tree::writeable_tree() {
    std::cout << "wtffffffff: " << Form("%s_vx", this->tag) << std::endl;
    this->tree->Branch(Form("%s_vx", this->tag), &this->vx);
    this->tree->Branch(Form("%s_vy", this->tag), &this->vy);
    this->tree->Branch(Form("%s_vz", this->tag), &this->vz);
    this->tree->Branch(Form("%s_vz_vpd", this->tag), &this->vpd_vz);
    
    this->tree->Branch(Form("%s_ep_east", this->tag), &this->ep_east);
    this->tree->Branch(Form("%s_ep_west", this->tag), &this->ep_west);

    this->tree->Branch(Form("%s_run_number", this->tag), &this->run_number);
    this->tree->Branch(Form("%s_num_triggers", this->tag), &this->num_triggers);
    this->tree->Branch(Form("%s_triggers", this->tag), this->triggers, Form("%s_triggers[%s_num_triggers]/L", this->tag, this->tag));

    this->tree->Branch(Form("%s_tofmult", this->tag), &this->tofmult);
    this->tree->Branch(Form("%s_tofmatch", this->tag), &this->tofmatch);
    this->tree->Branch(Form("%s_refmult3", this->tag), &this->refmult3);
    this->tree->Branch(Form("%s_centrality", this->tag), &this->centrality);
    this->tree->Branch(Form("%s_bbc_east_rate", this->tag), &this->bbc_east_rate);
    this->tree->Branch(Form("%s_bbc_west_rate", this->tag), &this->bbc_west_rate);

}

int Event_Tree::readable_tree() {}
