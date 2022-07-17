#include "tree_manager.h"

#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TBranch.h>

#include <map>
#include <string>
#include <iostream>

Tree_Manager::Tree_Manager() {};

Tree_Manager::Tree_Manager(TTree *_tree, std::string _tag) {
    this->tree = _tree;
    if (this->tree == nullptr) {
        throw std::runtime_error(Form("%s:%d tree is null", __FILE__, __LINE__));
    }
    this->tag = _tag.c_str();
};

Tree_Manager::~Tree_Manager() {};

void Tree_Manager::fill_tree() {
    if (this->tree == nullptr) {
        throw std::runtime_error(Form("%s:%d tree is null", __FILE__, __LINE__));
    }
    this->tree->Fill();
};
