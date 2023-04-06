#pragma once

class Zero_Crossing_Detector 
{
public:
    float xn = 0.0f;
    bool process(const float x);
};