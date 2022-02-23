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
#include "draw_histogram.h"

#include <iostream>

int main(int argc, char **argv) {
    TFile *file = TFile::Open("out.root");
    TTree *tree = file->Get<TTree>("jet_data");

    jet_tree_data datum;

    read_tree(tree, &datum);

    // Histograms
    histogram_package jet_spectrum;
    jet_spectrum.hist = new histogram_data;
    jet_spectrum.hist->hist = new TH1D("jet_spectra", "", 50, 0, 50);
    jet_spectrum.title = "Jet Spectra, Background Subtracted";
    jet_spectrum.x_title = "Momentum (GeV)";
    jet_spectrum.saveas = "plots/jet_spectra";


    histogram_data jet_eta;
    jet_eta.hist = new TH1D("jet_eta", "", 50, -1.5, 1.5);
    histogram_data jet_phi;
    jet_phi.hist = new TH1D("jet_phi", "", 50, 0, TMath::TwoPi());

    for (Long64_t n = 0; n < tree->GetEntries(); n++) {
        tree->GetEntry(n);

        // Fill histograms
        jet_spectrum.hist->hist->Fill(datum.jet_momentum_medium_subtracted);
        jet_eta.hist->Fill(datum.jet_eta);
        jet_phi.hist->Fill(datum.jet_phi);        
    }

    std::cerr << "made it here" << std::endl;

    draw_th1(&jet_spectrum);

}