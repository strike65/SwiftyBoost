#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "bs_complex.h"

// Spherical harmonics
bs_complex_d bs_spherical_harmonic(unsigned int n, int m, double theta, double phi);
bs_complex_f bs_spherical_harmonic_f(unsigned int n, int m, float theta, float phi);
bs_complex_l bs_spherical_harmonic_l(unsigned int n, int m, long double theta, long double phi);

#ifdef __cplusplus
}
#endif

