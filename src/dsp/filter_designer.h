#pragma once

struct FIR_Filter_Res {
    float* b;
    const int N;
    FIR_Filter_Res(const int _N): N(_N) 
    {
        b = new float[N];
    }
    ~FIR_Filter_Res() {
        delete [] b;
    }    
};

struct IIR_Filter_Res 
{
    float* a;
    float* b;
    const int N;
    IIR_Filter_Res(const int _N): N(_N) 
    {
        a = new float[N];
        b = new float[N];
    }
    ~IIR_Filter_Res() {
        delete [] a;
        delete [] b;
    }    
};

FIR_Filter_Res* create_fir_lpf(const float k, const int M);
FIR_Filter_Res* create_fir_hpf(const float k, const int M);
FIR_Filter_Res* create_fir_bpf(const float k1, const float k2, const int M);
IIR_Filter_Res* create_iir_single_pole_lpf(const float k);
IIR_Filter_Res* create_iir_notch_filter(const float k, const float r);
IIR_Filter_Res* create_iir_peak_filter(const float k, const float r);
