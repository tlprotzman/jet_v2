#include "event_tree.h"

#include "TROOT.h"
#include "TTree.h"

Event_Tree::Event_Tree(TTree *_tree, std::string _tag) {
    Tree_Manager(_tree, _tag);
    this->triggers = (long*)malloc(this->num_entries * sizeof(long));
}

Event_Tree::~Event_Tree() {
    free(this->triggers);
}

int Event_Tree::writeable_tree() {
    this->tree->Branch(Form("%s/vx", this->tag), &this->vx);
    this->tree->Branch(Form("%s/vy", this->tag), &this->vy);
    this->tree->Branch(Form("%s/vz", this->tag), &this->vz);
    this->tree->Branch(Form("%s/vz_vpd", this->tag), &this->vpd_vz);
    
    this->tree->Branch(Form("%s/ep_east", this->tag), &this->ep_east);
    this->tree->Branch(Form("%s/ep_west", this->tag), &this->ep_west);

    this->tree->Branch(Form("%s/run_number", this->tag), &this->run_number);
    this->tree->Branch(Form("%s/num_triggers", this->tag), &this->num_triggers);
    this->tree->Branch(Form("%s/triggers", this->tag), this->triggers, "triggers[num_triggers]/L");

    this->tree->Branch(Form("%s/tofmult", this->tag), &this->tofmult);
    this->tree->Branch(Form("%s/tofmatch", this->tag), &this->tofmatch);
    this->tree->Branch(Form("%s/refmult3", this->tag), &this->refmult3);
    this->tree->Branch(Form("%s/centrality", this->tag), &this->centrality);
    this->tree->Branch(Form("%s/bbc_east_rate", this->tag), &this->bbc_east_rate);
    this->tree->Branch(Form("%s/bbc_west_rate", this->tag), &this->bbc_west_rate);

}