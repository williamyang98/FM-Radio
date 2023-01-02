#pragma once

// Enable intrinsic code that can be compiled on target
#if defined(__AVX2__)
#define _DSP_AVX2
#endif

#if defined(__AVX2__) || defined(__SSSE3__)
#define _DSP_SSSE3
#endif

// On MSVC if __AVX2__ is defined then we have FMA
// On GCC __FMA__ is a given define
#if !defined(__FMA__) && defined(__AVX2__)
#define _DSP_FMA
#elif defined(__FMA__)
#define _DSP_FMA
#endif

#if defined(_DSP_AVX2)
#pragma message("Compiling DSP SIMD using AVX2 code")
#elif defined(_DSP_SSSE3)
#pragma message("Compiling DSP SIMD using SSSE3 code")
#else
#pragma message("Compiling DSP SIMD using scalar code")
#endif

#if defined(_DSP_FMA)
#pragma message("Compiling DSP SIMD with FMA code")
#endif