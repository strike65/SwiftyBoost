#pragma once

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

// Chebyshev
double bs_chebyshev_T(unsigned int n, double x);
double bs_chebyshev_U(unsigned int n, double x);
float bs_chebyshev_T_f(unsigned int n, float x);
float bs_chebyshev_U_f(unsigned int n, float x);
long double bs_chebyshev_T_l(unsigned int n, long double x);
long double bs_chebyshev_U_l(unsigned int n, long double x);

// Chebyshev series evaluation via Clenshaw (first kind)
double bs_chebyshev_clenshaw(const double *c, size_t count, double x);
float bs_chebyshev_clenshaw_f(const float *c, size_t count, float x);
long double bs_chebyshev_clenshaw_l(const long double *c, size_t count, long double x);

#ifdef __cplusplus
}
#endif

