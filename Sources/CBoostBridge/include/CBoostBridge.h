//
//  Created by VT on 11.10.25.
//  Copyright © 2025 Volker Thieme. All rights reserved.
//
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

// Boost.Math constants (long double) — only meaningful where Swift Float80 is available.
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

float bs_tgamma_f(float x);
float bs_lgamma_f(float x);
float bs_erf_f(float x);
float bs_erfc_f(float x);

long double bs_tgamma_l(long double x);
long double bs_lgamma_l(long double x);
long double bs_erf_l(long double x);
long double bs_erfc_l(long double x);

// Beta family
double bs_beta(double a, double b);
double bs_fullBeta(double a, double b, double x);
double bs_ibeta(double a, double b, double x);
double bs_ibetac(double a, double b, double x);
double bs_ibeta_inv(double a, double b, double p);
double bs_ibetac_inv(double a, double b, double p);
double bs_ibeta_inva(double b, double x, double p);
double bs_ibeta_invb(double a, double x, double p);

float bs_beta_f(float a, float b);
float bs_fullBeta_f(float a, float b, float x);
float bs_ibeta_f(float a, float b, float x);
float bs_ibetac_f(float a, float b, float x);
float bs_ibeta_inv_f(float a, float b, float p);
float bs_ibetac_inv_f(float a, float b, float p);
float bs_ibeta_inva_f(float b, float x, float p);
float bs_ibeta_invb_f(float a, float x, float p);

long double bs_beta_l(long double a, long double b);
long double bs_fullBeta_l(long double a, long double b, long double x);
long double bs_ibeta_l(long double a, long double b, long double x);
long double bs_ibetac_l(long double a, long double b, long double x);
long double bs_ibeta_inv_l(long double a, long double b, long double p);
long double bs_ibetac_inv_l(long double a, long double b, long double p);
long double bs_ibeta_inva_l(long double b, long double x, long double p);
long double bs_ibeta_invb_l(long double a, long double x, long double p);

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
double bs_expint_En(int n, double x);
double bs_expm1(double x);
double bs_log1p(double x);
double bs_log1pmx(double x);
double bs_powm1(double x, double y);
double bs_cbrt(double x);

float bs_expint_En_f(int n, float x);
float bs_expm1_f(float x);
float bs_log1p_f(float x);
float bs_log1pmx_f(float x);
float bs_powm1_f(float x, float y);
float bs_cbrt_f(float x);

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
double bs_ellint_3_complete(double k, double nu);
double bs_ellint_3(double k, double nu, double phi);

float bs_ellint_1_complete_f(float k);
float bs_ellint_1_f(float k, float phi);
float bs_ellint_2_complete_f(float k);
float bs_ellint_2_f(float k, float phi);
float bs_ellint_3_complete_f(float k, float nu);
float bs_ellint_3_f(float k, float nu, float phi);

long double bs_ellint_1_complete_l(long double k);
long double bs_ellint_1_l(long double k, long double phi);
long double bs_ellint_2_complete_l(long double k);
long double bs_ellint_2_l(long double k, long double phi);
long double bs_ellint_3_complete_l(long double k, long double nu);
long double bs_ellint_3_l(long double k, long double nu, long double phi);

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

#ifdef __cplusplus
}
#endif
