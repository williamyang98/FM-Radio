#pragma once

#include <complex>
#include "utility/span.h"

typedef struct fftwf_plan_s* fftwf_plan;

void CalculateFFT(
    tcb::span<const std::complex<float>> x,
    tcb::span<std::complex<float>> y);

void CalculateIFFT(
    tcb::span<const std::complex<float>> x,
    tcb::span<std::complex<float>> y);