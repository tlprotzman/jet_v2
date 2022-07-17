#include "histogram_data.h"

#include <TROOT.h>
#include <Rtypes.h>
#include <TH1.h>

histogram_data::histogram_data() {
    histogram_data(nullptr, "");
}

histogram_data::histogram_data(TH1 *hist, std::string name) {
    histogram_data(hist, name, kBlack);
}

histogram_data::histogram_data(TH1 *hist, std::string name, int color) {
    this->hist = hist;
    this->name = name;
    this->color = color;

    this->hist->SetLineColor(this->color);
}

histogram_data::~histogram_data() {
    return;
}

TH1 *histogram_data::get_hist() {
    return this->hist;
}

std::string histogram_data::get_name() {
    return this->name;
}

int histogram_data::get_color() {
    return this->color;
}