/*
Tristan Protzman
Lehigh University
March 10, 2022
tlprotzman@gmail.com
*/

#include "TROOT.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THStack.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TNamed.h"
#include "TStyle.h"
#include "TTree.h"

#include "setup.h"


int main(int argc, char **argv) {
    TFile *file = TFile::Open("out.root");
    TTree *tree = file->Get<TTree>("jet_data");
    jet_tree_data datum;
    read_tree(tree, &datum);

    TH1D *jet_relative_phi = new TH1D("jet_relative_phi", "Jet Spectra", 25, 0, TMath::PiOver2());
    
    // Populate relative angles
    double calculated_v2;
    for (Long64_t n = 0; n < tree->GetEntries(); n++) {
        tree->GetEvent(n);
        if (datum.centrality < 6 || datum.centrality > 9) {  // accept between 50-30% central
            continue;
        }

        for (int i = 0; i < datum.num_all_jets; i++) {
            double jet_phi = datum.all_jets_phi[i];
            double relative = jet_phi - datum.event_plane_full;
            if (relative < 0) {
                relative += TMath::TwoPi();
            }

            while (relative > TMath::PiOver2()) {
                relative -= TMath::PiOver2();
            }

            jet_relative_phi->Fill(relative);
            calculated_v2 += cos(2 * relative);
        }
    }
    calculated_v2 /= jet_relative_phi->GetEntries();
    std::cout << "Calculated V2: " << calculated_v2 << std::endl;

    // Calculate V2
    TF1 *v2_fit = new TF1("v2", "[0] * (1 + [1] * 2 * cos(2*x))", 0, TMath::PiOver2());
    v2_fit->SetParameters(1, 0);
    v2_fit->SetParNames("offset", "v2");
    jet_relative_phi->Fit("v2");



    TCanvas c("", "", 1000, 1000);

    jet_relative_phi->SetLineColor(kRed);



    std::cout << jet_relative_phi->GetFunction("v2")->GetParameter("v2") << std::endl;


    gPad->SetMargin(0.14, 0.9, 0.9, 0.9);
    jet_relative_phi->Draw("e");
    jet_relative_phi->SetXTitle("#Delta#phi=|#phi^{jet}-#Psi_{2}|");
    jet_relative_phi->SetYTitle("#frac{dN}{d#Delta#phi}");
    v2_fit->SetLineColor(kRed);
    jet_relative_phi->GetFunction("v2")->DrawF1(-1 * TMath::PiOver2(), TMath::PiOver2(), "C same");
    v2_fit->SetParameters(jet_relative_phi->GetFunction("v2")->GetParameter("offset"), calculated_v2);
    v2_fit->SetLineColor(kBlue);
    v2_fit->DrawF1(-1 * TMath::PiOver2(), TMath::PiOver2(), "same");
    TLatex formula = TLatex();
    formula.SetTextSize(0.02);
    // formula.DrawLatexNDC(0.35, 0.85, Form("Equation: Scale * (1 + 2v_{2} cos(2#phi))", jet_relative_phi->GetFunction("v2")->GetParameter("offset"), jet_relative_phi->GetFunction("v2")->GetParameter("v2")));
    formula.DrawLatexNDC(0.35, 0.85, Form(" red: %.1f(1 + 2 * %.5f cos(2#Delta#phi))", jet_relative_phi->GetFunction("v2")->GetParameter("offset"), jet_relative_phi->GetFunction("v2")->GetParameter("v2")));
    formula.DrawLatexNDC(0.35, 0.82, Form("blue: %.1f(1 + 2 * %.5f cos(2#Delta#phi))", jet_relative_phi->GetFunction("v2")->GetParameter("offset"), calculated_v2));
    formula.DrawLatexNDC(0.35, 0.79, "30\%-50\% central");
    
    formula.DrawLatexNDC(0.2, 0.2, "Red: Fitted to data");
    formula.DrawLatexNDC(0.2, 0.15, "Blue: v_{2}= #LTcos(2#Delta#phi)#GT");
    c.Print("plots/jet_relative_phi.png");
    c.Print("plots/jet_relative_phi.c");


    return 0;
}