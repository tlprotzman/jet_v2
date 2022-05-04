#include "tree_manager.h"

class Event_Tree: public Tree_Manager {
    public:
        // Vertex data
        double vx;
        double vy;
        double vz;
        double vpd_vz;
    
        // Event plane data
        double ep_east;
        double ep_west;

        // Triggers    
        int run_number;
        UInt_t num_triggers;
        long *triggers;

        // Detector Information
        unsigned short tofmult;
        unsigned short tofmatch;
        int refmult3;
        int centrality;
        float bbc_east_rate;
        float bbc_west_rate;     

        Event_Tree(TTree *_tree, std::string _tag);
        ~Event_Tree();

        int writeable_tree();
        int readable_tree();

        
};