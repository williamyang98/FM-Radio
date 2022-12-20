#pragma once

#include <complex>

class LocalOscillator 
{
public:
    float dt;
    float Ts;
    float f_center;
public:
    LocalOscillator();
    std::complex<float> Update();
};