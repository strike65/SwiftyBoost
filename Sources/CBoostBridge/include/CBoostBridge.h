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

// Boost.Math constants (as double)
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

double bs_const_euler(void);       // Euler–Mascheroni constant
double bs_const_catalan(void);     // Catalan's constant
double bs_const_zeta_three(void);  // Apery's constant ζ(3)
double bs_const_phi(void);         // Golden ratio φ

// Gamma / Error
double bs_tgamma(double x);
double bs_lgamma(double x);
double bs_erf(double x);
double bs_erfc(double x);

// Beta family
double bs_beta(double a, double b);
double bs_fullBeta(double a, double b, double x);
double bs_ibeta(double a, double b, double x);
double bs_ibetac(double a, double b, double x);
double bs_ibeta_inv(double a, double b, double p);
double bs_ibetac_inv(double a, double b, double p);
double bs_ibeta_inva(double b, double x, double p);
double bs_ibeta_invb(double a, double x, double p);

// Digamma / Polygamma / Zeta
double bs_digamma(double x);
double bs_trigamma(double x);
double bs_polygamma(int n, double x);
double bs_riemann_zeta(double x);

// Owen's T
double bs_owens_t(double h, double a);

// Exponential integrals and related
double bs_expint_En(int n, double x);
double bs_expm1(double x);
double bs_log1p(double x);
double bs_log1pmx(double x);
double bs_powm1(double x, double y);
double bs_cbrt(double x);

// Trig helpers
double bs_sin_pi(double x);
double bs_cos_pi(double x);

// Airy
double bs_airy_ai(double x);
double bs_airy_bi(double x);
double bs_airy_ai_prime(double x);
double bs_airy_bi_prime(double x);

// Bessel (cylindrical, real)
double bs_cyl_bessel_j(double v, double x);
double bs_cyl_neumann(double v, double x);
double bs_cyl_bessel_i(double v, double x);
double bs_cyl_bessel_k(double v, double x);

// Legendre (associated P only for integer n,m)
double bs_legendre_p(int n, double x);
double bs_assoc_legendre_p(int n, int m, double x);

// Elliptic integrals (Legendre forms)
double bs_ellint_1_complete(double k);
double bs_ellint_1(double k, double phi);
double bs_ellint_2_complete(double k);
double bs_ellint_2(double k, double phi);
double bs_ellint_3_complete(double k, double nu);
double bs_ellint_3(double k, double nu, double phi);

// Elliptic integrals (Carlson symmetric forms)
double bs_ellint_rc(double x, double y);
double bs_ellint_rf(double x, double y, double z);
double bs_ellint_rd(double x, double y, double z);
double bs_ellint_rj(double x, double y, double z, double p);
double bs_ellint_rg(double x, double y, double z);

// Lambert W
double bs_lambert_w0(double x);
double bs_lambert_wm1(double x);

#ifdef __cplusplus
}
#endif
