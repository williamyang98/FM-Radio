#pragma once

// Create an FIR filter with N taps
// b is a vector of length N
// k = Fc/(Fs/2)
void create_fir_lpf(float* b, const int N, const float k);
void create_fir_hpf(float* b, const int N, const float k);
void create_fir_bpf(float* b, const int N, const float k1, const float k2);

// Create a IIR single order buttworth LPF with 2 taps
// b, a are vectors of length 2 
// k = Fc/(Fs/2)
void create_iir_single_pole_lpf(float* b, float* a, const float k);

// Create a IIR second order peak/notch filter with 3 taps
// b, a are vectors of length 3
// k = Fc/(Fs/2)
void create_iir_notch_filter(float* b, float* a, const float k, const float r);
void create_iir_peak_filter(float* b, float* a, const float k, const float r);

// Create an FIR Hilbert filter with N taps
// b is a vector of length N
// For best results N should be odd
void create_fir_hilbert(float* b, const int N);

constexpr int TOTAL_TAPS_IIR_SINGLE_POLE_LPF = 2;
constexpr int TOTAL_TAPS_IIR_SECOND_ORDER_NOTCH_FILTER = 3;
constexpr int TOTAL_TAPS_IIR_SECOND_ORDER_PEAK_FILTER = 3;
