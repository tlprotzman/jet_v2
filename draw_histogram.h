#ifndef draw_histogram_h
#define draw_histogram_h

#include <stdlib.h>
#include "TH1.h"
#include "Rtypes.h"
#include <string>


typedef struct {
    TH1 *hist;
    std::string name="";
    EColor color=kBlack;
} histogram_data_s;

typedef struct {
    size_t num_histograms=1;
    histogram_data_s *hist;
    std::string title="";
    std::string x_title="";
    std::string y_title="";
    bool log_y=false;
    bool legend=false;
    std::string saveas;
    std::string format="png";
} histogram_package_s;

void draw_th1(histogram_package_s *histograms);
void draw_th1(histogram_data_s *histogram, std::string saveas);


#endif // draw_histogram_h