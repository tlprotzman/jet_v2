#include "particle_v2.h"

#include <TROOT.h>
#include <TH1.h>
#include <TH1D.h>
#include <TF1.h>

#include "fastjet/PseudoJet.hh"

#include <vector>

particle_v2 *calculate_particle_v2(std::vector<fastjet::PseudoJet> &tracks, fastjet::PseudoJet &leading_jet) {
    std::vector<fastjet::PseudoJet> background_tracks();
    for (fastjet::PseudoJet track : tracks) {
        if (true) {
            // background_tracks.push_back(track);
        }
    }
}
