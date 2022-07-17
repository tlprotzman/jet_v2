#ifndef HISTOGRAM_PACKAGE_H
#define HISTOGRAM_PACKAGE_H

#include "histogram_data.h"

#include <TROOT.h>
#include <Rtypes.h>
#include <TH1.h>
#include <TObject.h>

#include <string>

class histogram_package {
    public:
        // Constructor
        histogram_package();
        ~histogram_package();

        // Setters
        void set_title(std::string title);
        void set_x_title(std::string title);
        void set_y_title(std::string title);
        void set_save_location(std::string path);
        void set_format(std::string format);
        
        bool set_log_x(bool set);
        bool set_log_x() {return set_log_x(true);};
        bool set_log_y(bool set);
        bool set_log_y() {return set_log_y(true);};
        bool set_log_z(bool set);
        bool set_log_z() {return set_log_z(true);};
        bool set_legend(bool set);
        bool set_legend() {return set_legend(true);};
        bool set_legend(bool set, float x, float y);

        bool set_x_range(float x_min, float x_max);
        bool set_y_range(float y_min, float y_max);

        void add_histogram(histogram_data *hist);
        void add_drawable(TObject *drawable);

        void draw();

    // private:
        std::vector<histogram_data*> *histograms;
        std::vector<TObject*> *drawables;
        std::string title;
        std::string x_title;
        std::string y_title;
        
        bool log_x;
        bool log_y;
        bool log_z;
        
        bool legend;
        bool legend_man_pos;
        float legend_x, legend_y;

        bool man_x_range;
        float x_min, x_max;

        bool man_y_range;
        float y_min, y_max;

        std::string saveas;
        std::string format;
};

#endif // HISTOGRAM_PACKAGE_H
