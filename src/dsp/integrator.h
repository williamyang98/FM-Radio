#pragma once

template <typename T>
class Integrator_Block 
{
public:
    float KTs = 1.0f;
    T yn = 0.0f;
    T process(T x) {
        T y = KTs*x + yn; 
        yn = y;
        return y;
    }
};