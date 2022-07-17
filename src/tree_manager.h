#ifndef TREE_MANAGER_H
#define TREE_MANAGER_H

#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TBranch.h>

#include <map>
#include <string>

class Tree_Manager {
    protected:
        TTree *tree;
        const char *tag;

    public:
        const int num_entries = 1000;
        
        Tree_Manager();
        Tree_Manager(TTree *_tree, std::string _tag);
        ~Tree_Manager();

        virtual int writeable_tree() {};
        virtual int readable_tree() {};
        void fill_tree();
};

#endif // TREE_MANAGER_H