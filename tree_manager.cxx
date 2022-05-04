#include "tree_manager.h"

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"

#include <map>
#include <string>

Tree_Manager::Tree_Manager(TTree *_tree, std::string _tag) {
    this->tree = _tree;
    this->tag = _tag;
}

Tree_Manager::~Tree_Manager() {

}
