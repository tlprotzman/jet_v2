TCanvas *canvas = new TCanvas("", "", 1000, 1000);
    jet_momentum->Draw("hist");
    gPad->SetLogy();
    canvas->Print("plots/jet_momentum.png");
    jet_momentum_jet_median_subtracted->Draw("hist");
    gPad->SetLogy();
    canvas->Print("plots/jet_momentum_jet_subtracted.png");
    jet_momentum_grid_median_subtracted->Draw("hist");
    gPad->SetLogy();
    canvas->Print("plots/jet_momentum_grid_subtracted.png");
    track_constituent_momentum->Draw("hist");
    gPad->SetLogy();
    canvas->Print("plots/track_constituent_momentum.png");
    TCanvas *ep_canvas = new TCanvas("ep", "", 1000, 1000);
    for (uint32_t i = 0; i < 6; i+=2) {
        THStack *ep_stack = new THStack();
        TLegend *ep_legend = new TLegend();
        
        ep[i]->SetLineColor(kRed);
        ep_stack->Add(ep[i]);
        ep_legend->AddEntry(ep[i], "East");
        
        ep[i+1]->SetLineColor(kBlue);
        ep_stack->Add(ep[i+1]);
        ep_legend->AddEntry(ep[i+1], "West");
        
        ep_stack->Draw("nostack");
        ep_stack->GetXaxis()->SetTitle("Phi");
        ep_stack->GetYaxis()->SetTitle("Counts");
        ep_stack->SetTitle(ep_hist_names[i] + 5);
        ep_legend->Draw();
        ep_canvas->Print(Form("plots/%s.png", ep_hist_names[i] + 5));
        delete ep_stack;
        delete ep_legend;
    }