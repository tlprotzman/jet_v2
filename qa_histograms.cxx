/*
Tristan Protzman
Lehigh University
February 22, 2022
tlprotzman@gmail.com

Produces (hopefully) pretty QA plots from the jet v2 analysis
*/

#include "TROOT.h"
#include "TTree.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH1D.h"
#include "THStack.h"
#include "TLegend.h"
#include "TFile.h"
#include "TMath.h"
#include "TStyle.h"

#include "setup.h"
#include "histogram_data.h"
#include "histogram_package.h"

#include <iostream>

int main(int argc, char **argv) {
    TFile *file = TFile::Open("out.root");
    TTree *tree = file->Get<TTree>("jet_data");

    jet_tree_data datum;

    read_tree(tree, &datum);

    // Histograms
    histogram_package jet_pt_spectrum;
    jet_pt_spectrum.set_save_location("plots/nocut_pt_spectrum");
    jet_pt_spectrum.set_log_y(true);
    jet_pt_spectrum.set_title("Jet Momentum Spectra");
    jet_pt_spectrum.set_x_title("Pt");
    jet_pt_spectrum.set_y_title("Count");

    jet_pt_spectrum.add_histogram(new histogram_data(file->Get<TH1>("jet_momentum"), "Jet Mometum", kBlue));
    jet_pt_spectrum.add_histogram(new histogram_data(file->Get<TH1>("jet_subtracted_momentum"), "Jet Mometum, Background Subtracted", kRed));

    jet_pt_spectrum.set_legend(true, 0.48, 0.7);
    jet_pt_spectrum.set_x_range(-15, 25);
    
    jet_pt_spectrum.draw();
}