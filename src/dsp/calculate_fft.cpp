#include "calculate_fft.h"

#define _USE_MATH_DEFINES
#include <cmath>
#include <fftw3.h>
#include <mutex>
#include <unordered_map>

struct Key 
{ 
    size_t block_size; 
    bool is_inverse; 
    bool operator==(const Key& other) const {
        return 
            (block_size == other.block_size) &&
            (is_inverse == other.is_inverse);
    }
};

struct KeyHasher 
{
    std::size_t operator()(const Key& k) const {
        const size_t shift = sizeof(size_t)*8 - 1;
        return (k.block_size | ((size_t)k.is_inverse << shift));
    }
};

static auto fft_plans = std::unordered_map<Key, fftwf_plan, KeyHasher>();
static auto mutex_fft_plans = std::mutex();

static fftwf_plan GetPlan(const size_t block_size, const bool is_inverse) {
    auto lock = std::scoped_lock(mutex_fft_plans);
    auto key = Key{ block_size, is_inverse };
    auto res = fft_plans.find(key);
    if (res == fft_plans.end()) {
        auto type = is_inverse ? FFTW_BACKWARD : FFTW_FORWARD;
        auto plan = fftwf_plan_dft_1d((int)block_size, NULL, NULL, type, FFTW_ESTIMATE);
        res = fft_plans.insert({ key, plan }).first;
    }
    return res->second;
}

void CalculateFFT(
    tcb::span<const std::complex<float>> x,
    tcb::span<std::complex<float>> y)
{
    const size_t N = x.size();
    auto plan = GetPlan(N, false);
    fftwf_execute_dft(plan, (fftwf_complex*)x.data(), (fftwf_complex*)y.data());
}

void CalculateIFFT(
    tcb::span<const std::complex<float>> x,
    tcb::span<std::complex<float>> y)
{
    const size_t N = x.size();
    auto plan = GetPlan(N, true);
    fftwf_execute_dft(plan, (fftwf_complex*)x.data(), (fftwf_complex*)y.data());
}