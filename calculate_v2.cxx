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
#include "TGraphErrors.h"
#include "TNamed.h"
#include "TStyle.h"
#include "TTree.h"

#include "setup.h"
#include "histogram_data.h"
#include "histogram_package.h"


int main(int argc, char **argv) {
    // ROOT::EnableImplicitMT();
    TFile *file = TFile::Open("out.root");
    TTree *tree = file->Get<TTree>("jet_data");
    jet_tree_data datum;
    read_tree(tree, &datum);

    size_t num_pt_bins = 8;
    double pt_bin_upper_bound[8] = {5, 10, 15, 20, 25, 30, 35, 40};
    std::vector<TH1D*> jet_relative_phi;
    std::vector<TH1D*> jet_relative_phi_unfolded;
    std::vector<double> calculated_v2;

    for (size_t i = 0; i < num_pt_bins + 1; i++) {
        jet_relative_phi.push_back(new TH1D(Form("jet_relative_phi_%i", i), "Jet Spectra", 25, 0, TMath::PiOver2()));
        jet_relative_phi_unfolded.push_back(new TH1D(Form("jet_relative_phi_unfolded_%i", i), "Jet Spectra", 100, 0, TMath::TwoPi()));
        calculated_v2.push_back(0);
    }
    
    // Populate relative angles
    for (Long64_t n = 0; n < tree->GetEntries(); n++) {
        // if (n > 1000000) {
        //     break;
        // }
        if (n % 100000 == 0) {
            std::cout << "Processed " << n << " events" << std::endl;
        }
        tree->GetEvent(n);
        if (datum.centrality < 6 || datum.centrality > 9) {  // accept between 50-30% central
            continue;
        }


        for (int i = 0; i < datum.num_all_jets; i++) {
            // if (datum.all_jets_pt[i] < 10) {
            //     continue;
            // }
            if (datum.all_jets_z[i] > 0.95) {
                continue;
            }
            double jet_phi = datum.all_jets_phi[i];
            double relative = abs(jet_phi - datum.event_plane_full);
            double other_relative = jet_phi - datum.event_plane_full;
            if (other_relative < 0) {
                other_relative += TMath::TwoPi();
            }
            int pt_bin = 0;
            while (pt_bin < num_pt_bins && datum.all_jets_pt[i] > pt_bin_upper_bound[pt_bin]) {pt_bin++;}   // Populate the appropriate momentum bin
            jet_relative_phi_unfolded[pt_bin]->Fill(other_relative);
            // if (relative < 0) {
            //     relative += TMath::TwoPi();
            // }

            while (relative > TMath::PiOver2()) {
                relative -= TMath::PiOver2();
            }

            jet_relative_phi[pt_bin]->Fill(relative);
            calculated_v2[pt_bin] += cos(2 * relative);
        }
    }
    for (size_t i = 0; i < num_pt_bins + 1; i++) {
        calculated_v2[i] /= jet_relative_phi[i]->GetEntries();
        std::cout << "Bin " << i << " calculated V2: " << calculated_v2[i] << " from " << jet_relative_phi[i]->GetEntries() << " events" << std::endl;
    }

    // Calculate V2
    double *x = (double*)malloc(num_pt_bins + 1 * sizeof(double));
    double *y  = (double*)malloc(num_pt_bins * sizeof(double));
    double *ex = (double*)malloc(num_pt_bins * sizeof(double));
    double *ey = (double*)malloc(num_pt_bins * sizeof(double));
    std::vector<TF1*> v2_fits;
    for (size_t i = 0; i < num_pt_bins + 1; i++) {
        if (jet_relative_phi[i]->GetEntries() == 0) {
            continue;
        }
        TF1 *v2_fit = new TF1(Form("v2_bin_%d", i), "[0] * (1 + [1] * 2 * cos(2*x))", 0, TMath::PiOver2());
        v2_fit->SetParameters(1, 0);
        v2_fit->SetParNames("offset", "v2");
        jet_relative_phi[i]->Fit(Form("v2_bin_%d", i));
        x[i] = i * 5;
        ex[i] = 1;
        y[i] = v2_fit->GetParameter("v2");
        ey[i] = v2_fit->GetParError(v2_fit->GetParNumber("v2"));
        v2_fits.push_back(v2_fit);
    }

    TCanvas *c = new TCanvas("", "", 1000, 1000);
    c->SetMargin(0.15, 0.1, 0.1, 0.1);
    TGraphErrors *v2_momentum_dependence = new TGraphErrors(num_pt_bins, x, y, ex, ey);
    v2_momentum_dependence->GetXaxis()->SetTitle("P_{t} Bin");
    v2_momentum_dependence->GetYaxis()->SetTitle("v_{2}");
    v2_momentum_dependence->SetTitle("Momentum Dependence of v_{2}");
    v2_momentum_dependence->SetMarkerSize(2.5);
    v2_momentum_dependence->SetMarkerStyle(1);
    v2_momentum_dependence->SetLineWidth(2);
    v2_momentum_dependence->Draw("AP");
    c->SaveAs("plots/v2_momentum_dependence.png");
    c->SaveAs("plots/v2_momentum_dependence.c");
    return 0;

    



    /*
    histogram_package *v2_plots = new histogram_package();
    histogram_data **

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

    jet_relative_phi_unfolded->Draw("e");
    v2_full_range->SetLineColor(kRed);
    jet_relative_phi_unfolded->GetFunction("v2_full")->DrawF1(0, TMath::TwoPi(), "C same");
    v2_full_range->SetParameters(jet_relative_phi_unfolded->GetFunction("v2_full")->GetParameter("offset"), calculated_v2);
    v2_full_range->SetLineColor(kBlue);
    v2_full_range->DrawF1(0, TMath::TwoPi(), "same");
    // formula.DrawLatexNDC(0.35, 0.85, Form("Equation: Scale * (1 + 2v_{2} cos(2#phi))", jet_relative_phi_unfolded->GetFunction("v2")->GetParameter("offset"), jet_relative_phi_unfolded->GetFunction("v2")->GetParameter("v2")));
    formula.DrawLatexNDC(0.35, 0.85, Form(" red: %.1f(1 + 2 * %.5f cos(2#Delta#phi))", jet_relative_phi_unfolded->GetFunction("v2_full")->GetParameter("offset"), jet_relative_phi->GetFunction("v2")->GetParameter("v2")));
    formula.DrawLatexNDC(0.35, 0.82, Form("blue: %.1f(1 + 2 * %.5f cos(2#Delta#phi))", jet_relative_phi_unfolded->GetFunction("v2_full")->GetParameter("offset"), calculated_v2));
    formula.DrawLatexNDC(0.35, 0.79, "30\%-50\% central");
    
    formula.DrawLatexNDC(0.2, 0.2, "Red: Fitted to data");
    formula.DrawLatexNDC(0.2, 0.15, "Blue: v_{2}= #LTcos(2#Delta#phi)#GT");

    c.Print("plots/jet_relative_phi_unfolded.png");
    c.Print("plots/jet_relative_phi_unfolded.c");

    return 0;
    */
}