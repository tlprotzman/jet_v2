#include "histogram_package.h"
#include "histogram_data.h"

#include <TROOT.h>
#include <Rtypes.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <THStack.h>
#include <TObject.h>
#include <TStyle.h>

#include <string>
#include <iostream>

histogram_package::histogram_package() {
    this->histograms = new std::vector<histogram_data*>;
    this->drawables = new std::vector<TObject*>;
    this->title = "";
    this->x_title = "";
    this->y_title = "";

    this->log_x = false;
    this->log_y = false;
    this->log_z = false;
    this->legend = false;
    this->legend_man_pos = false;
    this->man_x_range = false;
    this->man_y_range = false;

    this->saveas = "";
    this->format = "png";
}

histogram_package::~histogram_package() {
    delete this->histograms;
}

void histogram_package::set_title(std::string title) {
    this->title = title;
}

void histogram_package::set_x_title(std::string title) {
    this->x_title = title;
}

void histogram_package::set_y_title(std::string title) {
    this->y_title = title;
}

void histogram_package::set_save_location(std::string path) {
    this->saveas = path;
}

void histogram_package::set_format(std::string format) {
    this->format = format;
}

bool histogram_package::set_log_x(bool set) {
    this->log_x = set;
    return this->log_x;
}

bool histogram_package::set_log_y(bool set) {
    this->log_y = set;
    return this->log_y;
}

bool histogram_package::set_log_z(bool set) {
    this->log_z = set;
    return this->log_z;
}

bool histogram_package::set_legend(bool set) {
    this->legend = set;
    return this->legend;
}

bool histogram_package::set_legend(bool set, float x, float y) {
    this->legend = set;
    this->legend_man_pos = true;
    this->legend_x = x;
    this->legend_y = y;
    return this->legend;
}

bool histogram_package::set_x_range(float x_min, float x_max) {
    this->man_x_range = true;
    this->x_min = x_min;
    this->x_max = x_max;
}

bool histogram_package::set_y_range(float y_min, float y_max) {
    this->man_y_range = true;
    this->y_min = y_min;
    this->y_max = y_max;
}

void histogram_package::add_histogram(histogram_data *hist) {
    this->histograms->push_back(hist);
}

void histogram_package::add_drawable(TObject *drawable) {
    this->drawables->push_back(drawable);
}

void histogram_package::draw() {
    TCanvas *c = new TCanvas("", "", 1000, 1000);
    THStack *stack = new THStack();
    TLegend *legend = new TLegend();

    c->SetMargin(0.155, 0.1, 0.1, 0.1);

    for (auto hist : *this->histograms) {
        stack->Add(hist->get_hist());
        legend->AddEntry(hist->get_hist(), hist->get_name().c_str());
    }
    stack->Draw("nostack");

    stack->SetTitle(this->title.c_str());
    stack->GetXaxis()->SetTitle(this->x_title.c_str());
    stack->GetYaxis()->SetTitle(this->y_title.c_str());

    gPad->SetLogx(this->log_x);
    gPad->SetLogy(this->log_y);
    gPad->SetLogz(this->log_z);

    if (this->man_x_range) {
        stack->GetXaxis()->SetRangeUser(this->x_min, this->x_max);
        stack->Draw("nostack");
    }

    if (this->man_y_range) {
        std::cout << "running" << std::endl;
        // stack->GetYaxis()->SetRange(this->y_min, this->y_max);
        stack->SetMinimum(this->y_min);
        stack->SetMaximum(this->y_max);
        stack->Draw("nostack");
    }

    if (this->legend) {
        legend->Draw();
        if (this->legend_man_pos) {
            gPad->Update();
            legend->SetX1NDC(this->legend_x);
            legend->SetX2NDC(this->legend_x + 0.4);
            legend->SetY1NDC(this->legend_y);
            legend->SetY2NDC(this->legend_y + 0.1);
            gPad->Modified();
            std::cout << "ran" << std::endl;
        }
    }

    // Draw remaining drawables
    for (auto drawable : *this->drawables) {
        drawable->Draw();
    }

    if (this->saveas != "") {
        c->SaveAs(Form("%s.%s", this->saveas.c_str(), this->format.c_str()));
        c->SaveAs(Form("%s.%s", this->saveas.c_str(), "c"));
    }
}
