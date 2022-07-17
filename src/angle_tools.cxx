#include "angle_tools.h"

#include <TMath.h>

double smallest_angle(double a, double b) {
    double relative_angle = a - b;
    if (relative_angle < 0) {
        relative_angle += TMath::TwoPi();
    }
    if (relative_angle > TMath::Pi()) { // Two angles can never be more than pi apart
        relative_angle -= TMath::Pi();
    }
    return relative_angle;
}