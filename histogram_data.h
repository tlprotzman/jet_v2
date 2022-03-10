#ifndef HISTOGRAM_DATA_H
#define HISTOGRAM_DATA_H


#include "TROOT.h"
#include "Rtypes.h"
#include "TH1.h"

class histogram_data {
    public:
        histogram_data();
        histogram_data(TH1 *hist, std::string name);
        histogram_data(TH1 *hist, std::string name, EColor color);
        ~histogram_data();

        TH1 *get_hist();
        std::string get_name();
        EColor get_color();

    private:
        TH1 *hist;
        std::string name;
        EColor color;
};

#endif // HISTOGRAM_DATA_H