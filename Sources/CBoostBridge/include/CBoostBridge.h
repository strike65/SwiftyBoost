//
//  Created by VT on 11.10.25.
//  Copyright © 2025 Volker Thieme.
//  Permission is hereby granted, free of charge, to any person obtaining a copy
//  of this software and associated documentation files (the "Software"), to deal
//  in the Software without restriction, including without limitation the rights
//  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//  copies of the Software, and to permit persons to whom the Software is
//  furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included in
//  all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
//  THE SOFTWARE.
//  

#pragma once
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

// Boost.Math constants (double)
double bs_const_e(void);
double bs_const_pi(void);
double bs_const_two_pi(void);
double bs_const_half_pi(void);
double bs_const_third_pi(void);
double bs_const_two_thirds_pi(void);
double bs_const_quarter_pi(void);
double bs_const_three_quarters_pi(void);
double bs_const_sixth_pi(void);
double bs_const_pi_sqr(void);

double bs_const_root_two(void);
double bs_const_root_three(void);
double bs_const_root_pi(void);

double bs_const_one_div_pi(void);
double bs_const_one_div_two_pi(void);
double bs_const_one_div_root_pi(void);
double bs_const_two_div_pi(void);
double bs_const_two_div_root_pi(void);

double bs_const_ln_two(void);
double bs_const_ln_ten(void);
double bs_const_ln_ln_two(void);

double bs_const_euler(void);
double bs_const_catalan(void);
double bs_const_zeta_three(void);
double bs_const_phi(void);

// Boost.Math constants (float)
float bs_const_e_f(void);
float bs_const_pi_f(void);
float bs_const_two_pi_f(void);
float bs_const_half_pi_f(void);
float bs_const_third_pi_f(void);
float bs_const_two_thirds_pi_f(void);
float bs_const_quarter_pi_f(void);
float bs_const_three_quarters_pi_f(void);
float bs_const_sixth_pi_f(void);
float bs_const_pi_sqr_f(void);

float bs_const_root_two_f(void);
float bs_const_root_three_f(void);
float bs_const_root_pi_f(void);

float bs_const_one_div_pi_f(void);
float bs_const_one_div_two_pi_f(void);
float bs_const_one_div_root_pi_f(void);
float bs_const_two_div_pi_f(void);
float bs_const_two_div_root_pi_f(void);

float bs_const_ln_two_f(void);
float bs_const_ln_ten_f(void);
float bs_const_ln_ln_two_f(void);

float bs_const_euler_f(void);
float bs_const_catalan_f(void);
float bs_const_zeta_three_f(void);
float bs_const_phi_f(void);

// Boost.Math constants (long double)
long double bs_const_e_l(void);
long double bs_const_pi_l(void);
long double bs_const_two_pi_l(void);
long double bs_const_half_pi_l(void);
long double bs_const_third_pi_l(void);
long double bs_const_two_thirds_pi_l(void);
long double bs_const_quarter_pi_l(void);
long double bs_const_three_quarters_pi_l(void);
long double bs_const_sixth_pi_l(void);
long double bs_const_pi_sqr_l(void);

long double bs_const_root_two_l(void);
long double bs_const_root_three_l(void);
long double bs_const_root_pi_l(void);

long double bs_const_one_div_pi_l(void);
long double bs_const_one_div_two_pi_l(void);
long double bs_const_one_div_root_pi_l(void);
long double bs_const_two_div_pi_l(void);
long double bs_const_two_div_root_pi_l(void);

long double bs_const_ln_two_l(void);
long double bs_const_ln_ten_l(void);
long double bs_const_ln_ln_two_l(void);

long double bs_const_euler_l(void);
long double bs_const_catalan_l(void);
long double bs_const_zeta_three_l(void);
long double bs_const_phi_l(void);

// Gamma / Error
double bs_tgamma(double x);
double bs_lgamma(double x);
double bs_erf(double x);
double bs_erfc(double x);
// Γ(a) / Γ(b)
double      bs_tgamma_ratio(double a, double b);
float       bs_tgamma_ratio_f(float a, float b);
long double bs_tgamma_ratio_l(long double a, long double b);

// Γ(a) / Γ(a + delta)
double      bs_tgamma_delta_ratio(double a, double delta);
float       bs_tgamma_delta_ratio_f(float a, float delta);
long double bs_tgamma_delta_ratio_l(long double a, long double delta);

float bs_tgamma_f(float x);
float bs_lgamma_f(float x);
float bs_erf_f(float x);
float bs_erfc_f(float x);

long double bs_tgamma_l(long double x);
long double bs_lgamma_l(long double x);
long double bs_erf_l(long double x);
long double bs_erfc_l(long double x);

// Incomplete gamma (lower/upper, regularized, and inverses)
double bs_tgamma_lower(double a, double x);
double bs_tgamma_upper(double a, double x);
double bs_gamma_p(double a, double x);
double bs_gamma_q(double a, double x);
double bs_gamma_p_inv(double a, double p);
double bs_gamma_q_inv(double a, double q);

// Derivatives of regularized incomplete gamma (w.r.t. x)
double bs_gamma_p_derivative(double a, double x);

float bs_tgamma_lower_f(float a, float x);
float bs_tgamma_upper_f(float a, float x);
float bs_gamma_p_f(float a, float x);
float bs_gamma_q_f(float a, float x);
float bs_gamma_p_inv_f(float a, float p);
float bs_gamma_q_inv_f(float a, float q);

// Derivatives (float)
float bs_gamma_p_derivative_f(float a, float x);

long double bs_tgamma_lower_l(long double a, long double x);
long double bs_tgamma_upper_l(long double a, long double x);
long double bs_gamma_p_l(long double a, long double x);
long double bs_gamma_q_l(long double a, long double x);
long double bs_gamma_p_inv_l(long double a, long double p);
long double bs_gamma_q_inv_l(long double a, long double q);

// Derivatives (long double)
long double bs_gamma_p_derivative_l(long double a, long double x);

// Beta family
double bs_beta(double a, double b);
double bs_fullBeta(double a, double b, double x);
double bs_ibeta(double a, double b, double x);
double bs_ibetac(double a, double b, double x);
double bs_ibeta_inv(double a, double b, double p);
double bs_ibetac_inv(double a, double b, double p);
double bs_ibeta_inva(double b, double x, double p);
double bs_ibeta_invb(double a, double x, double p);
// Derivative of regularized incomplete Beta I_x(a, b) w.r.t. x
double bs_ibeta_derivative(double a, double b, double x);

float bs_beta_f(float a, float b);
float bs_fullBeta_f(float a, float b, float x);
float bs_ibeta_f(float a, float b, float x);
float bs_ibetac_f(float a, float b, float x);
float bs_ibeta_inv_f(float a, float b, float p);
float bs_ibetac_inv_f(float a, float b, float p);
float bs_ibeta_inva_f(float b, float x, float p);
float bs_ibeta_invb_f(float a, float x, float p);
// Derivative (float)
float bs_ibeta_derivative_f(float a, float b, float x);

long double bs_beta_l(long double a, long double b);
long double bs_fullBeta_l(long double a, long double b, long double x);
long double bs_ibeta_l(long double a, long double b, long double x);
long double bs_ibetac_l(long double a, long double b, long double x);
long double bs_ibeta_inv_l(long double a, long double b, long double p);
long double bs_ibetac_inv_l(long double a, long double b, long double p);
long double bs_ibeta_inva_l(long double b, long double x, long double p);
long double bs_ibeta_invb_l(long double a, long double x, long double p);
// Derivative (long double)
long double bs_ibeta_derivative_l(long double a, long double b, long double x);

// Digamma / Polygamma / Zeta
double bs_digamma(double x);
double bs_trigamma(double x);
double bs_polygamma(int n, double x);
double bs_riemann_zeta(double x);

float bs_digamma_f(float x);
float bs_trigamma_f(float x);
float bs_polygamma_f(int n, float x);
float bs_riemann_zeta_f(float x);

long double bs_digamma_l(long double x);
long double bs_trigamma_l(long double x);
long double bs_polygamma_l(int n, long double x);
long double bs_riemann_zeta_l(long double x);

// Owen's T
double bs_owens_t(double h, double a);
float bs_owens_t_f(float h, float a);
long double bs_owens_t_l(long double h, long double a);

// Exponential integrals and related
double bs_expint_Ei(double x);
double bs_expint_En(int n, double x);
double bs_expm1(double x);
double bs_log1p(double x);
double bs_log1pmx(double x);
double bs_powm1(double x, double y);
double bs_cbrt(double x);

float bs_expint_Ei_f(float x);
float bs_expint_En_f(int n, float x);
float bs_expm1_f(float x);
float bs_log1p_f(float x);
float bs_log1pmx_f(float x);
float bs_powm1_f(float x, float y);
float bs_cbrt_f(float x);

long double bs_expint_Ei_l(long double x);
long double bs_expint_En_l(int n, long double x);
long double bs_expm1_l(long double x);
long double bs_log1p_l(long double x);
long double bs_log1pmx_l(long double x);
long double bs_powm1_l(long double x, long double y);
long double bs_cbrt_l(long double x);

// Trig helpers
double bs_sin_pi(double x);
double bs_cos_pi(double x);

float bs_sin_pi_f(float x);
float bs_cos_pi_f(float x);

long double bs_sin_pi_l(long double x);
long double bs_cos_pi_l(long double x);

// Airy
double bs_airy_ai(double x);
double bs_airy_bi(double x);
double bs_airy_ai_prime(double x);
double bs_airy_bi_prime(double x);

float bs_airy_ai_f(float x);
float bs_airy_bi_f(float x);
float bs_airy_ai_prime_f(float x);
float bs_airy_bi_prime_f(float x);

long double bs_airy_ai_l(long double x);
long double bs_airy_bi_l(long double x);
long double bs_airy_ai_prime_l(long double x);
long double bs_airy_bi_prime_l(long double x);

// Bessel (cylindrical, real)
double bs_cyl_bessel_j(double v, double x);
double bs_cyl_neumann(double v, double x);
double bs_cyl_bessel_i(double v, double x);
double bs_cyl_bessel_k(double v, double x);

float bs_cyl_bessel_j_f(float v, float x);
float bs_cyl_neumann_f(float v, float x);
float bs_cyl_bessel_i_f(float v, float x);
float bs_cyl_bessel_k_f(float v, float x);

long double bs_cyl_bessel_j_l(long double v, long double x);
long double bs_cyl_neumann_l(long double v, long double x);
long double bs_cyl_bessel_i_l(long double v, long double x);
long double bs_cyl_bessel_k_l(long double v, long double x);

// Legendre
double bs_legendre_p(int n, double x);
double bs_assoc_legendre_p(int n, int m, double x);

float bs_legendre_p_f(int n, float x);
float bs_assoc_legendre_p_f(int n, int m, float x);

long double bs_legendre_p_l(int n, long double x);
long double bs_assoc_legendre_p_l(int n, int m, long double x);

// Elliptic integrals (Legendre forms)
double bs_ellint_1_complete(double k);
double bs_ellint_1(double k, double phi);
double bs_ellint_2_complete(double k);
double bs_ellint_2(double k, double phi);
double bs_ellint_3(double k, double nu, double phi); // incomplete elliptic integral of the third kind ∏(k, nu, phi)
double bs_ellint_3_complete(double k, double nu);

float bs_ellint_1_complete_f(float k);
float bs_ellint_1_f(float k, float phi);
float bs_ellint_2_complete_f(float k);
float bs_ellint_2_f(float k, float phi);
float bs_ellint_3_f(float k, float nu, float phi); // incomplete elliptic integral of the third kind ∏(k, nu, phi)
float bs_ellint_3_complete_f(float k, float nu);

long double bs_ellint_1_complete_l(long double k);
long double bs_ellint_1_l(long double k, long double phi);
long double bs_ellint_2_complete_l(long double k);
long double bs_ellint_2_l(long double k, long double phi);
long double bs_ellint_3_l(long double k, long double nu,long double phi); // incomplete elliptic integral of the third kind ∏(k, nu, phi)
long double bs_ellint_3_complete_l(long double k, long double nu);

// Elliptic integrals (Carlson symmetric forms)
double bs_ellint_rc(double x, double y);
double bs_ellint_rf(double x, double y, double z);
double bs_ellint_rd(double x, double y, double z);
double bs_ellint_rj(double x, double y, double z, double p);
double bs_ellint_rg(double x, double y, double z);

float bs_ellint_rc_f(float x, float y);
float bs_ellint_rf_f(float x, float y, float z);
float bs_ellint_rd_f(float x, float y, float z);
float bs_ellint_rj_f(float x, float y, float z, float p);
float bs_ellint_rg_f(float x, float y, float z);

long double bs_ellint_rc_l(long double x, long double y);
long double bs_ellint_rf_l(long double x, long double y, long double z);
long double bs_ellint_rd_l(long double x, long double y, long double z);
long double bs_ellint_rj_l(long double x, long double y, long double z, long double p);
long double bs_ellint_rg_l(long double x, long double y, long double z);

// Lambert W
double bs_lambert_w0(double x);
double bs_lambert_wm1(double x);

float bs_lambert_w0_f(float x);
float bs_lambert_wm1_f(float x);

long double bs_lambert_w0_l(long double x);
long double bs_lambert_wm1_l(long double x);

// Hypergeometric functions (Gauss/confluent/general)
double bs_hypergeometric_1F0(double a, double z);
float  bs_hypergeometric_1F0_f(float a, float z);
long double bs_hypergeometric_1F0_l(long double a, long double z);

double bs_hypergeometric_0F1(double b, double z);
float  bs_hypergeometric_0F1_f(float b, float z);
long double bs_hypergeometric_0F1_l(long double b, long double z);

double bs_hypergeometric_2F0(double a, double b, double z);
float  bs_hypergeometric_2F0_f(float a, float b, float z);
long double bs_hypergeometric_2F0_l(long double a, long double b, long double z);

double bs_hypergeometric_1F1(double a, double b, double z);
float  bs_hypergeometric_1F1_f(float a, float b, float z);
long double bs_hypergeometric_1F1_l(long double a, long double b, long double z);

// General pFq with arrays of parameters
double bs_hypergeometric_pFq(const double* a, size_t p, const double* b, size_t q, double z);
float  bs_hypergeometric_pFq_f(const float* a, size_t p, const float* b, size_t q, float z);
long double bs_hypergeometric_pFq_l(const long double* a, size_t p, const long double* b, size_t q, long double z);

// Bernoulli Numbers
double bs_bernoulli_b2n(const int n);
float  bs_bernoulli_b2n_f(const int n);
long double bs_bernoulli_b2n_l(const int n);

// Tangent Numbers (scalar)
double bs_tangent_t2n(const int n);
float  bs_tangent_t2n_f(const int n);
long double bs_tangent_t2n_l(const int n);

// Tangent Numbers (bulk sequence)
// Fills out[count] with T_{2(start_index)}, T_{2(start_index+1)}, ..., count values.
void bs_tangent_t2n_seq(int start_index, unsigned int count, double* out);
void bs_tangent_t2n_seq_f(int start_index, unsigned int count, float* out);
// long double; // keep C header valid for long double type in prototypes below
void bs_tangent_t2n_seq_l(int start_index, unsigned int count, long double* out);

// prime numbers up to 10000th
unsigned int bs_prime(unsigned int n);

// Fibonacci numbers
// Returns the nth Fibonacci number as unsigned long long (F(0)=0, F(1)=1).
unsigned long long bs_fibonacci_ull(unsigned long long n);

// Factorials
double bs_factorial(unsigned int i);
float bs_factorial_f(unsigned int i);
long double bs_factorial_l(unsigned int i);

// Pochhammer
double bs_rising_factorial(double x, unsigned int i);
float bs_rising_factorial_f(float x, unsigned int i);
long double bs_rising_factorial_l(long double x, unsigned int i);

// binomial coefficients
double bs_binomial_coefficient( unsigned int n, unsigned int k);
float bs_binomial_coefficient_f( unsigned int n, unsigned int k);
long double bs_binomial_coefficient_l( unsigned int n, unsigned int k);

// double factorial
double bs_double_factorial(unsigned int i);
float bs_double_factorial_f(unsigned int i);
long double bs_double_factorial_l(unsigned int i);

#ifdef __cplusplus
}
#endif

