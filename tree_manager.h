#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"

#include <map>
#include <string>

class Tree_Manager {
    protected:
        TTree *tree;
        std::string tag;

    public:
        const int num_entries = 10;
        
        Tree_Manager();
        Tree_Manager(TTree *_tree, std::string _tag);
        ~Tree_Manager();

        virtual int writeable_tree();
        virtual int readable_tree();
};