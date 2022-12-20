#include "zero_crossing_detector.h"

bool Zero_Crossing_Detector::process(float x) {
    // crossing occurs if the signals are on opposite sides of the x-axis
    bool is_crossed = (x*xn) < 0.0f;
    xn = x;
    return is_crossed;
}