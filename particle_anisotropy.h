#ifndef PARTICLE_V2_H
#define PARTICLE_V2_H

#include "fastjet/PseudoJet.hh"

#include <TH1.h>

#include <vector>

struct particle_anisotropy {
    double rho;
    double v2;
    double v3;
};

void calculate_v2(std::vector<fastjet::PseudoJet> &tracks, fastjet::PseudoJet &leading_jet, double event_plane, particle_anisotropy &flow_parameters, int centrality);
void print_particle_distribution(TH1 *hist);

#endif // PARTICLE_V2_H