// Throws particles with an assumed v2/v3 and tries to recover v2/v3,
// both with a hole in phi and without to test effects of missing TPC sector

#include <TROOT.h>

#include <TH1.h>
#include <TH1D.h>
#include <TF1.h>
#include <TTree.h>
#include <TMath.h>
#include <TRandom.h>
#include <TRandom3.h>

#include <iostream>

#include "tree_manager.h"
#include "jet_tree.h"
#include "event_tree.h"

int main(int argc, char **argv) {
    int n = 100000;
    const double v2 = 0.05;
    const double v3 = 0.00;
    const double multiplicity = 300;
    const double percent_missing = 1;

    TFile *outfile = new TFile("simout.root", "RECREATE");
    TTree *tree = new TTree("jet_data", "Jet Data");
    tree->SetDirectory(outfile);
    Event_Tree *event_tree = new Event_Tree(tree, "event");
    event_tree->writeable_tree();    
    Jet_Tree *jet_tree = new Jet_Tree(tree, "all");
    jet_tree->writeable_tree();

    TRandom *rng = new TRandom3();
    rng->SetSeed(4); // Selected by roll of fair die

    TF1 *distribution = new TF1("azimuthal", "[0] + 2 * ([1] * cos(2 * (x - [3])) + [2] * cos(3 * (x - [3])))", 0, TMath::TwoPi());
    distribution->SetParName(0, "Amplitude");
    distribution->SetParName(1, "v2");
    distribution->SetParName(2, "v3");
    distribution->SetParName(3, "event_plane");

    distribution->SetParameter("Amplitude", 1);
    distribution->SetParameter("v2", v2);
    distribution->SetParameter("v3", v3);


    for (int i = 0; i < n; i++) {
        if (i % 1000 == 0) {
            std::cout << "Processed " << i << " events" << std::endl;
        }
        double ep = rng->Rndm() * TMath::TwoPi();
        distribution->SetParameter("event_plane", ep);
        event_tree->ep_east = ep + 0.2;
        event_tree->ep_west = ep + 0.2;
        int j = 0;
        while (j < multiplicity) {
            double angle = rng->Rndm() * TMath::TwoPi();
            double prob = rng->Rndm() * 1.5;
            double is_missing = rng->Rndm();
            if (prob < distribution->Eval(angle)) {
                if (angle > 4 && angle < 5 && is_missing < percent_missing) {
                    continue;
                }
                jet_tree->jet_phi[j] = angle;
                j++;
            }
        }
        jet_tree->num_jets = j;
        event_tree->fill_tree();
    }

    std::cout << "Saving..." << std::endl;
    tree->Write();
    outfile->Close();
    return 0;
}
