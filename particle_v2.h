#ifndef PARTICLE_V2_H
#define PARTICLE_V2_H

#include <TROOT.h>

#include "fastjet/PseudoJet.hh"

#include <vector>

typedef struct particle_v2 {
    double rho;
    double v2;
    double v3;
};

particle_v2 *calculate_particle_v2(std::vector<fastjet::PseudoJet> &tracks, fastjet::PseudoJet &leading_jet);






#endif //PARTICLE_V2_H