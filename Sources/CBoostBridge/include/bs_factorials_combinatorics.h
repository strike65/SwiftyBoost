#pragma once

#ifdef __cplusplus
extern "C" {
#endif

// Factorials
double bs_factorial(unsigned int i);
float bs_factorial_f(unsigned int i);
long double bs_factorial_l(unsigned int i);

// Pochhammer (rising factorial)
double bs_rising_factorial(double x, unsigned int i);
float bs_rising_factorial_f(float x, unsigned int i);
long double bs_rising_factorial_l(long double x, unsigned int i);

// Binomial coefficients
double bs_binomial_coefficient(unsigned int n, unsigned int k);
float bs_binomial_coefficient_f(unsigned int n, unsigned int k);
long double bs_binomial_coefficient_l(unsigned int n, unsigned int k);

// Double factorial
double bs_double_factorial(unsigned int i);
float bs_double_factorial_f(unsigned int i);
long double bs_double_factorial_l(unsigned int i);

#ifdef __cplusplus
}
#endif

