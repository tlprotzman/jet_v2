/*
Tristan Protzman
Lehigh University
February 23, 2022
tlprotzman@gmail.com

Let's see if I can make a reasonably generic plotting macro that I'm happy with...
*/

#include "draw_histogram.h"

#include "TH1.h"
#include "THStack.h"
#include "TLegend.h"
#include "TCanvas.h"


void draw_th1(histogram_package *histograms) {
    TCanvas *c = new TCanvas("", "", 1000, 1000);
    THStack *stack = new THStack();
    TLegend *legend = new TLegend();

    for (size_t i = 0; i < histograms->num_histograms; i++) {
        histogram_data hist = histograms->hist[i];
        hist.hist->SetLineColor(hist.color);
        stack->Add(hist.hist);
        legend->AddEntry(hist.hist, hist.name.c_str());
    }

    // stack->SetTitle(histograms->title.c_str());
    // stack->GetXaxis()->SetTitle(histograms->x_title.c_str());
    // stack->GetYaxis()->SetTitle(histograms->y_title.c_str());

    stack->Draw("nostack");
    if (histograms->legend) {
        legend->Draw();
    }
    c->SaveAs(Form("%s.%s", histograms->saveas.c_str(), histograms->format.c_str()));
    c->SaveAs(Form("%s.%s", histograms->saveas.c_str(), "c"));
}

void draw_th1(histogram_data *histogram, std::string saveas) {
    histogram_package package;
    package.hist = histogram;
    package.saveas = saveas;
    draw_th1(&package);
}