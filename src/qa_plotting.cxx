#include <TROOT.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TGraph.h>

#include <TCanvas.h>
#include <TLegend.h>
#include <THStack.h>
#include <TStyle.h>
#include <TLatex.h>

void qa_plotting() {
    gStyle->SetOptStat(0);
    TFile *infile = new TFile("testout.root");
    TLatex *text = new TLatex();
    text->SetTextSize(0.037);
    text->SetTextFont(42);

    TH1 *jet_spectra = infile->Get<TH1>("jet_pt_spectra");
    TH1 *jet_spectra_no_match = infile->Get<TH1>("jet_pt_spectra_nomatch");
    TH1 *hc_spectra = infile->Get<TH1>("hc_pt_spectra");

    THStack *pt_spectra = new THStack();
    TLegend *legend = new TLegend(0.60, 0.52, 0.895, 0.70);
    legend->SetTextFont(42);
    legend->SetTextSize(0.038);
    legend->SetBorderSize(0);

    
    jet_spectra->SetLineColor(kRed);
    jet_spectra->SetFillColorAlpha(kRed, 0.6);
    jet_spectra->SetFillStyle(3003);
    jet_spectra_no_match->SetLineColor(kGreen+2);
    jet_spectra_no_match->SetFillColorAlpha(kGreen+1, 0.7);
    jet_spectra_no_match->SetFillStyle(3003);
    hc_spectra->SetLineColor(kBlue);
    hc_spectra->SetFillColorAlpha(kBlue+1, 0.4);
    hc_spectra->SetFillStyle(3002);
    hc_spectra->GetXaxis()->SetRangeUser(10, 50);

    TH1 *jet_spectra_points = (TH1*)jet_spectra->Clone("jet_spectra_points");
    TH1 *jet_spectra_no_match_points = (TH1*)jet_spectra_no_match->Clone("jet_spectra_no_match_points");
    TH1 *hc_spectra_points = (TH1*)hc_spectra->Clone("hc_spectra_points");

    jet_spectra->SetLineWidth(0);
    jet_spectra_no_match->SetLineWidth(0);
    hc_spectra->SetLineWidth(0);

    jet_spectra_points->SetLineWidth(3);
    jet_spectra_no_match_points->SetLineWidth(3);
    hc_spectra_points->SetLineWidth(3);

    jet_spectra_no_match->GetXaxis()->SetRangeUser(10, 50);
    double scale_factor = 1.0 / jet_spectra_no_match->Integral();
    jet_spectra_no_match->Scale(scale_factor);
    jet_spectra->Scale(scale_factor);
    hc_spectra->Scale(scale_factor);
    
    jet_spectra_points->Scale(scale_factor);
    jet_spectra_no_match_points->Scale(scale_factor);
    hc_spectra_points->Scale(scale_factor);
    
    pt_spectra->Add(jet_spectra_no_match);
    pt_spectra->Add(jet_spectra);
    pt_spectra->Add(hc_spectra);

    

    legend->AddEntry(jet_spectra_no_match_points, "Jets");
    legend->AddEntry(hc_spectra_points, "Hard Cores");
    legend->AddEntry(jet_spectra_points, "Matched Jets");

    TCanvas *c = new TCanvas("c", "c", 1000, 1000);
    c->SetLeftMargin(0.125);
    pt_spectra->Draw("nostack hist");
    pt_spectra->GetXaxis()->SetTitle("p_{T, jet}^{reco, ch}");
    pt_spectra->GetXaxis()->SetTitleSize(0.038);
    pt_spectra->GetXaxis()->SetTitleOffset(1.15);
    pt_spectra->GetYaxis()->SetTitle("Normalized Counts");
    pt_spectra->GetYaxis()->SetTitleSize(0.038);
    pt_spectra->GetXaxis()->SetRangeUser(10, 50);
    legend->Draw();
    jet_spectra_no_match_points->Draw("e same");
    hc_spectra_points->Draw("e same");
    jet_spectra_points->Draw("e same");
    text->SetTextFont(62);
    text->SetTextFont(42);
    text->DrawLatexNDC(0.79, 0.86, "STAR");
    // text->DrawLatexNDC(0.635, 0.71, "Preliminary");
    text->DrawLatexNDC(0.512, 0.82, "Jet Spectra, Uncorrected");
    text->DrawLatexNDC(0.43, 0.78, "#sqrt{s_{NN}}=200 GeV Ru+Ru & Zr+Zr");
    text->DrawLatexNDC(0.682, 0.74, "Anti-k_{T} R=0.2");
    text->DrawLatexNDC(0.766, 0.70, "|#eta| < 0.8");
    c->SetLogy();

    c->SaveAs("plots/qa_pt_spectra.png");
    c->SaveAs("plots/qa_pt_spectra.pdf");

    TH3 *matching_performance = infile->Get<TH3>("matching_performance");
    TH1 *deta_dphi = matching_performance->Project3D("xy");
    deta_dphi->SetTitle("Match Jet Distance");
    TCanvas *c2 = new TCanvas("c2", "c2", 1000, 1000);
    c2->SetMargin(0.125, 0.125, 0.125, 0.125);
    deta_dphi->Draw("colz");
    c2->SetLogz();
    c2->SaveAs("plots/qa_matched_distance.png");
    c2->SaveAs("plots/qa_matched_distance.pdf");

    TCanvas *c3 = new TCanvas("c3", "c3", 1000, 1000);
    TH1 *dpt = matching_performance->Project3D("z");
    dpt->SetTitle("Match Jet p_{T} Difference");
    dpt->GetXaxis()->SetTitle("Jet p_{T} - HC p_{T}");
    dpt->GetXaxis()->SetRangeUser(-15, 15);
    dpt->Draw("hist");
    // c3->SetLogy();
    c3->SaveAs("plots/qa_matched_pt.png");
    c3->SaveAs("plots/qa_matched_pt.pdf");



    TH2 *ue_subtraction = infile->Get<TH2>("ue_subtraction");
    ue_subtraction->SetTitle("Underlying Event");
    TCanvas *c4 = new TCanvas("c4", "", 1000, 1000);
    c4->SetMargin(0.125, 0.125, 0.125, 0.125);
    ue_subtraction->GetYaxis()->SetRangeUser(0, 8);
    ue_subtraction->Draw("colz");
    // c4->SetLogz();
    c4->SaveAs("plots/qa_ue_subtracted.png");
    c4->SaveAs("plots/qa_ue_subtracted.pdf");

    TH1 *relative_angle;
    c4->SetLeftMargin(0.18);
    c4->SetRightMargin(0.05);
    c4->SetTopMargin(0.05);
    std::vector<std::string> prefixes = {"hc", "all", "hc_only"};
    for (std::string prefix : prefixes) {
        int slice = 1;
        while ((relative_angle = infile->Get<TH1>(Form("relative_angle_%s_%i", prefix.c_str(), slice))) != nullptr) {
            double scale_factor = 1.0 / relative_angle->Integral();
            relative_angle->Scale(scale_factor);
            relative_angle->SetLineWidth(2);
            TF1 *func = relative_angle->GetFunction(relative_angle->GetListOfFunctions()->First()->GetName());
            // std::cout << "offset was " << func->GetParameter("offset") << std::endl;
            func->SetParameter("offset", func->GetParameter("offset") * scale_factor);
            // std::cout << "offset is " << func->GetParameter("offset") << std::endl;
            TLegend *l = new TLegend(0.2, 0.15, 0.5, 0.3);
            l->SetTextSize(0.038);
            l->SetLineWidth(0);
            l->AddEntry(relative_angle, "Jet Yield");
            TH1 *boxes_are_weird = new TH1D("linesareweird", "", 1, 1, 1);
            boxes_are_weird->SetLineWidth(2);
            boxes_are_weird->SetLineColor(kRed);
            l->AddEntry(boxes_are_weird, "#it{A}(1 + 2v_{2}^{obs}cos(2#Delta#phi)) Fit");

            relative_angle->Draw("e");
            l->Draw();

            text->DrawLatexNDC(0.66, 0.90, "STAR Preliminary");
            // text->DrawLatexNDC(0.845, 0.90, "STAR");
            text->DrawLatexNDC(0.47, 0.86, "#sqrt{s_{NN}}=200 GeV Ru+Ru & Zr+Zr");
            text->DrawLatexNDC(0.729, 0.82, "Anti-k_{T} R=0.2");
            text->DrawLatexNDC(0.813, 0.78, "|#eta| < 0.8");
            // text->DrawLatexNDC(0.702, 0.74, relative_angle->GetTitle() + 11);
            text->DrawLatexNDC(0.55, 0.74, "10 < p_{T, jet}^{reco, ch} < 12.5 GeV/c");
            // text->DrawLatexNDC(0.632, 0.70, "Momentum Uncorrected");
            text->DrawLatexNDC(0.62, 0.69, "Statistical Error Only");
            // text->DrawLatexNDC(0.655, 0.66, Form("v_{2}^{ch, obs} = %.3f #pm %.3f", func->GetParameter("v2"), func->GetParError(func->GetParNumber("v2"))));

            relative_angle->GetXaxis()->SetTitleOffset(1.2);
            relative_angle->GetXaxis()->SetTitle("#Delta#phi = #Psi_{2} - #phi_{jet}");
            relative_angle->GetYaxis()->SetTitle("#frac{1}{N} #frac{d^{2}N}{dp_{T} d#Delta#phi}");
            relative_angle->GetYaxis()->SetTitleOffset(2.3);
            relative_angle->SetTitle("");
            c4->SaveAs(Form("plots/qa_slice_%s_%i.pdf", prefix.c_str(), slice));
            c4->SaveAs(Form("plots/qa_slice_%s_%i.png", prefix.c_str(), slice));
            delete l;
            slice++;
        }
    }


    infile->Close();

    TCanvas *c5 = new TCanvas("c5", "c5", 1600, 1000);
    c5->SetLeftMargin(0.1);
    c5->SetTopMargin(0.05);
    c5->SetRightMargin(0.05);
    TFile *infile_2 = new TFile("/data/star/jet_v2/8_10_R02/histograms_r02");
    
    TH1 *east_uncorrected = infile_2->Get<TH1>("east_uncorrected");
    TH1 *east_phi_corrected = infile_2->Get<TH1>("east_phi_corrected");
    TH1 *east_phi_psi_corrected = infile_2->Get<TH1>("east_phi_psi_corrected");

    // east_uncorrected->SetBinContent(1, east_uncorrected->GetBinContent(2) * 0.995);
    // east_phi_corrected->SetBinContent(1, east_phi_corrected->GetBinContent(2) * 0.995);
    // east_phi_psi_corrected->SetBinContent(1, east_phi_psi_corrected->GetBinContent(2) * 0.995);

    THStack *ep_flattening = new THStack();
    TLegend *ep_legend = new TLegend(0.12, 0.12, 0.55, 0.3);
    ep_legend->SetTextSize(0.03);
    ep_legend->SetLineWidth(0);

    east_uncorrected->SetLineColor(kRed);
    east_uncorrected->Scale(1.0 / east_uncorrected->Integral());
    east_uncorrected->SetLineWidth(2);
    ep_flattening->Add(east_uncorrected);
    ep_legend->AddEntry(east_uncorrected, "Uncorrected");

    east_phi_corrected->SetLineColor(kBlue);
    east_phi_corrected->Scale(1.0 / east_phi_corrected->Integral());
    east_phi_corrected->SetLineWidth(2);
    ep_flattening->Add(east_phi_corrected);
    ep_legend->AddEntry(east_phi_corrected, "Phi Weighted");

    east_phi_psi_corrected->SetLineColor(kBlack);
    east_phi_psi_corrected->Scale(1.0 / east_phi_psi_corrected->Integral());
    east_phi_psi_corrected->SetLineWidth(2);
    ep_flattening->Add(east_phi_psi_corrected);
    ep_legend->AddEntry(east_phi_psi_corrected, "Phi Weighted, Psi Shifted");

    ep_flattening->Draw("nostack");
    // ep_flattening->GetYaxis()->SetRangeUser(0.03, 0.04);
    ep_flattening->SetMinimum(0.03);
    ep_flattening->SetMaximum(0.036);
    ep_flattening->GetXaxis()->SetTitle("#phi");
    ep_flattening->GetYaxis()->SetTitle("Event Fraction");
    ep_flattening->Draw("nostack");

    ep_legend->Draw();

    text->DrawLatexNDC(0.13, 0.40, "STAR");
    text->DrawLatexNDC(0.13, 0.36, "#sqrt{s_{NN}}=200 GeV Ru+Ru & Zr+Zr");
    text->DrawLatexNDC(0.13, 0.32, "East EPD #Psi_{2}");

    c5->SaveAs("plots/ep_flattening.png");
    c5->SaveAs("plots/ep_flattening.pdf");

    TCanvas *c6 = new TCanvas("c6", "c6", 1600, 1000);
    c6->SetGrid();
    c6->SetRightMargin(0.05);
    c6->SetTopMargin(0.05);
    TH2 *ep_resolution = infile_2->Get<TH2>("ep_resolution");
    TGraph *ep_reso_graph = new TGraph();
    // TH1 *dummy_hist = new TH1I("dummy", "", 0, 16, 1);
    // dummy_hist->Draw("AXIS");

    for (int i = 1; i <= 16; i++) {
        ep_resolution->GetYaxis()->SetRange(i, i);
        TH1 *reso_projection = ep_resolution->ProjectionX(Form("ep_reso_%i", i));
        double reso = 0;
        int n = 0;
        for (int j = 1; j <= reso_projection->GetXaxis()->GetNbins(); j++) {
            // std::cout << reso_projection->GetBinCenter(j) << "\t" << reso_projection->GetBinContent(j) << std::endl;
            reso += cos(2 * reso_projection->GetBinCenter(j)) * reso_projection->GetBinContent(j);
            n += reso_projection->GetBinContent(j);
        }
        reso /= n;
        reso = sqrt(2 * reso);
        std::cout << "Resolution bin " << i << ": " << reso << std::endl;
        ep_reso_graph->SetPoint(i-1, i - 0.5, reso);

        // reso_projection->Draw("hist");
        // c5->SaveAs(Form("plots/ep_resolution_new_%i.pdf", i));
        // c5->SaveAs(Form("plots/ep_resolution_new_%i.png", i));
    }

    ep_reso_graph->SetMarkerStyle(kFullCircle);
    ep_reso_graph->SetMarkerColor(kBlack);
    ep_reso_graph->SetMarkerSize(2);
    ep_reso_graph->SetLineWidth(0);

    ep_reso_graph->Draw("");
    // dummy_hist->GetYaxis()->SetRangeUser(0, 0.5);
    ep_reso_graph->GetYaxis()->SetRangeUser(0, 0.45);
    ep_reso_graph->GetXaxis()->SetRangeUser(0, 15.9);
    ep_reso_graph->GetXaxis()->SetLabelOffset(100);

    for (uint i = 0; i <= 8; i++) {
        // text->SetTextAngle(90);
        text->DrawText((i*2)-0.2, -0.02, Form("%i%%", i * 10));
    }

    ep_reso_graph->GetYaxis()->SetTitle("#Psi_{2} Resolution");
    ep_reso_graph->GetYaxis()->SetTitleSize(0.04);
    ep_reso_graph->GetXaxis()->SetTitle("Centrality");
    ep_reso_graph->GetXaxis()->SetTitleSize(0.04);
    ep_reso_graph->GetXaxis()->SetTitleOffset(1.15);


    text->SetTextSize(0.04);
    text->DrawLatexNDC(0.88, 0.90, "STAR");
    text->DrawLatexNDC(0.75, 0.86, "EPD #Psi_{2} Resolution");
    text->DrawLatexNDC(0.64, 0.81, "#sqrt{s_{NN}}=200 GeV Ru+Ru & Zr+Zr");

    c6->SaveAs("plots/ep_resolution_new.pdf");   
    c6->SaveAs("plots/ep_resolution_new.png");


}