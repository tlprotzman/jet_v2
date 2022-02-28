#ifndef DICT_H
#define DICT_H

#pragma link C++ class std::vector<single_jet_data> +;


typedef struct {
    double pt;
    double e;
    double eta, phi;
    double subtracted_momentum;
    size_t num_constituents;
    double z;
} single_jet_data;

#endif // DICT_H