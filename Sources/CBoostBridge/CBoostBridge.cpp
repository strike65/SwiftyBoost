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

#include "CBoostBridge.h"

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>
#include <boost/math/special_functions/polygamma.hpp>
#include <boost/math/special_functions/zeta.hpp>
#include <boost/math/special_functions/owens_t.hpp>
#include <boost/math/special_functions/expint.hpp>
#include <boost/math/special_functions/expm1.hpp>
#include <boost/math/special_functions/log1p.hpp>
#include <boost/math/special_functions/powm1.hpp>
#include <boost/math/special_functions/cbrt.hpp>
#include <boost/math/special_functions/sinc.hpp>
#include <boost/math/special_functions/sinhc.hpp>
#include <boost/math/special_functions/sin_pi.hpp>
#include <boost/math/special_functions/cos_pi.hpp>
#include <boost/math/special_functions/airy.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include <boost/math/special_functions/ellint_3.hpp>
#include <boost/math/special_functions/ellint_rc.hpp>
#include <boost/math/special_functions/ellint_rf.hpp>
#include <boost/math/special_functions/ellint_rd.hpp>
#include <boost/math/special_functions/ellint_rj.hpp>
#include <boost/math/special_functions/ellint_rg.hpp>
#include <boost/math/special_functions/lambert_w.hpp>
#include <boost/math/constants/constants.hpp>

extern "C" {

// Boost constants (double)
double bs_const_e(void)                 { return boost::math::constants::e<double>(); }
double bs_const_pi(void)                { return boost::math::constants::pi<double>(); }
double bs_const_two_pi(void)            { return boost::math::constants::two_pi<double>(); }
double bs_const_half_pi(void)           { return boost::math::constants::half_pi<double>(); }
double bs_const_third_pi(void)          { return boost::math::constants::third_pi<double>(); }
double bs_const_two_thirds_pi(void)     { return boost::math::constants::two_thirds_pi<double>(); }
double bs_const_quarter_pi(void)        { return boost::math::constants::quarter_pi<double>(); }
double bs_const_three_quarters_pi(void) { return boost::math::constants::three_quarters_pi<double>(); }
double bs_const_sixth_pi(void)          { return boost::math::constants::sixth_pi<double>(); }
double bs_const_pi_sqr(void)            { return boost::math::constants::pi_sqr<double>(); }

double bs_const_root_two(void)          { return boost::math::constants::root_two<double>(); }
double bs_const_root_three(void)        { return boost::math::constants::root_three<double>(); }
double bs_const_root_pi(void)           { return boost::math::constants::root_pi<double>(); }

double bs_const_one_div_pi(void)        { return boost::math::constants::one_div_pi<double>(); }
double bs_const_one_div_two_pi(void)    { return boost::math::constants::one_div_two_pi<double>(); }
double bs_const_one_div_root_pi(void)   { return boost::math::constants::one_div_root_pi<double>(); }
double bs_const_two_div_pi(void)        { return boost::math::constants::two_div_pi<double>(); }
double bs_const_two_div_root_pi(void)   { return boost::math::constants::two_div_root_pi<double>(); }

double bs_const_ln_two(void)            { return boost::math::constants::ln_two<double>(); }
double bs_const_ln_ten(void)            { return boost::math::constants::ln_ten<double>(); }
double bs_const_ln_ln_two(void)         { return boost::math::constants::ln_ln_two<double>(); }

double bs_const_euler(void)             { return boost::math::constants::euler<double>(); }
double bs_const_catalan(void)           { return boost::math::constants::catalan<double>(); }
double bs_const_zeta_three(void)        { return boost::math::constants::zeta_three<double>(); }
double bs_const_phi(void)               { return boost::math::constants::phi<double>(); }

// Boost constants (float)
float bs_const_e_f(void)                 { return boost::math::constants::e<float>(); }
float bs_const_pi_f(void)                { return boost::math::constants::pi<float>(); }
float bs_const_two_pi_f(void)            { return boost::math::constants::two_pi<float>(); }
float bs_const_half_pi_f(void)           { return boost::math::constants::half_pi<float>(); }
float bs_const_third_pi_f(void)          { return boost::math::constants::third_pi<float>(); }
float bs_const_two_thirds_pi_f(void)     { return boost::math::constants::two_thirds_pi<float>(); }
float bs_const_quarter_pi_f(void)        { return boost::math::constants::quarter_pi<float>(); }
float bs_const_three_quarters_pi_f(void) { return boost::math::constants::three_quarters_pi<float>(); }
float bs_const_sixth_pi_f(void)          { return boost::math::constants::sixth_pi<float>(); }
float bs_const_pi_sqr_f(void)            { return boost::math::constants::pi_sqr<float>(); }

float bs_const_root_two_f(void)          { return boost::math::constants::root_two<float>(); }
float bs_const_root_three_f(void)        { return boost::math::constants::root_three<float>(); }
float bs_const_root_pi_f(void)           { return boost::math::constants::root_pi<float>(); }

float bs_const_one_div_pi_f(void)        { return boost::math::constants::one_div_pi<float>(); }
float bs_const_one_div_two_pi_f(void)    { return boost::math::constants::one_div_two_pi<float>(); }
float bs_const_one_div_root_pi_f(void)   { return boost::math::constants::one_div_root_pi<float>(); }
float bs_const_two_div_pi_f(void)        { return boost::math::constants::two_div_pi<float>(); }
float bs_const_two_div_root_pi_f(void)   { return boost::math::constants::two_div_root_pi<float>(); }

float bs_const_ln_two_f(void)            { return boost::math::constants::ln_two<float>(); }
float bs_const_ln_ten_f(void)            { return boost::math::constants::ln_ten<float>(); }
float bs_const_ln_ln_two_f(void)         { return boost::math::constants::ln_ln_two<float>(); }

float bs_const_euler_f(void)             { return boost::math::constants::euler<float>(); }
float bs_const_catalan_f(void)           { return boost::math::constants::catalan<float>(); }
float bs_const_zeta_three_f(void)        { return boost::math::constants::zeta_three<float>(); }
float bs_const_phi_f(void)               { return boost::math::constants::phi<float>(); }

// Boost constants (long double)
long double bs_const_e_l(void)                 { return boost::math::constants::e<long double>(); }
long double bs_const_pi_l(void)                { return boost::math::constants::pi<long double>(); }
long double bs_const_two_pi_l(void)            { return boost::math::constants::two_pi<long double>(); }
long double bs_const_half_pi_l(void)           { return boost::math::constants::half_pi<long double>(); }
long double bs_const_third_pi_l(void)          { return boost::math::constants::third_pi<long double>(); }
long double bs_const_two_thirds_pi_l(void)     { return boost::math::constants::two_thirds_pi<long double>(); }
long double bs_const_quarter_pi_l(void)        { return boost::math::constants::quarter_pi<long double>(); }
long double bs_const_three_quarters_pi_l(void) { return boost::math::constants::three_quarters_pi<long double>(); }
long double bs_const_sixth_pi_l(void)          { return boost::math::constants::sixth_pi<long double>(); }
long double bs_const_pi_sqr_l(void)            { return boost::math::constants::pi_sqr<long double>(); }

long double bs_const_root_two_l(void)          { return boost::math::constants::root_two<long double>(); }
long double bs_const_root_three_l(void)        { return boost::math::constants::root_three<long double>(); }
long double bs_const_root_pi_l(void)           { return boost::math::constants::root_pi<long double>(); }

long double bs_const_one_div_pi_l(void)        { return boost::math::constants::one_div_pi<long double>(); }
long double bs_const_one_div_two_pi_l(void)    { return boost::math::constants::one_div_two_pi<long double>(); }
long double bs_const_one_div_root_pi_l(void)   { return boost::math::constants::one_div_root_pi<long double>(); }
long double bs_const_two_div_pi_l(void)        { return boost::math::constants::two_div_pi<long double>(); }
long double bs_const_two_div_root_pi_l(void)   { return boost::math::constants::two_div_root_pi<long double>(); }

long double bs_const_ln_two_l(void)            { return boost::math::constants::ln_two<long double>(); }
long double bs_const_ln_ten_l(void)            { return boost::math::constants::ln_ten<long double>(); }
long double bs_const_ln_ln_two_l(void)         { return boost::math::constants::ln_ln_two<long double>(); }

long double bs_const_euler_l(void)             { return boost::math::constants::euler<long double>(); }
long double bs_const_catalan_l(void)           { return boost::math::constants::catalan<long double>(); }
long double bs_const_zeta_three_l(void)        { return boost::math::constants::zeta_three<long double>(); }
long double bs_const_phi_l(void)               { return boost::math::constants::phi<long double>(); }

// Gamma / Error
double bs_tgamma(double x)              { return boost::math::tgamma(x); }
double bs_lgamma(double x)              { return boost::math::lgamma(x); }
double bs_erf(double x)                 { return boost::math::erf(x); }
double bs_erfc(double x)                { return boost::math::erfc(x); }

float bs_tgamma_f(float x)              { return boost::math::tgamma(x); }
float bs_lgamma_f(float x)              { return boost::math::lgamma(x); }
float bs_erf_f(float x)                 { return boost::math::erf(x); }
float bs_erfc_f(float x)                { return boost::math::erfc(x); }

long double bs_tgamma_l(long double x)  { return boost::math::tgamma(x); }
long double bs_lgamma_l(long double x)  { return boost::math::lgamma(x); }
long double bs_erf_l(long double x)     { return boost::math::erf(x); }
long double bs_erfc_l(long double x)    { return boost::math::erfc(x); }

// Beta family
double bs_beta(double a, double b)                     { return boost::math::beta(a, b); }
double bs_fullBeta(double a, double b, double x)       { return boost::math::beta(a, b, x); }
double bs_ibeta(double a, double b, double x)          { return boost::math::ibeta(a, b, x); }
double bs_ibetac(double a, double b, double x)         { return boost::math::ibetac(a, b, x); }
double bs_ibeta_inv(double a, double b, double p)      { return boost::math::ibeta_inv(a, b, p); }
double bs_ibetac_inv(double a, double b, double p)     { return boost::math::ibetac_inv(a, b, p); }
double bs_ibeta_inva(double b, double x, double p)     { return boost::math::ibeta_inva(b, x, p); }
double bs_ibeta_invb(double a, double x, double p)     { return boost::math::ibeta_invb(a, x, p); }

float bs_beta_f(float a, float b)                      { return boost::math::beta(a, b); }
float bs_fullBeta_f(float a, float b, float x)         { return boost::math::beta(a, b, x); }
float bs_ibeta_f(float a, float b, float x)            { return boost::math::ibeta(a, b, x); }
float bs_ibetac_f(float a, float b, float x)           { return boost::math::ibetac(a, b, x); }
float bs_ibeta_inv_f(float a, float b, float p)        { return boost::math::ibeta_inv(a, b, p); }
float bs_ibetac_inv_f(float a, float b, float p)       { return boost::math::ibetac_inv(a, b, p); }
float bs_ibeta_inva_f(float b, float x, float p)       { return boost::math::ibeta_inva(b, x, p); }
float bs_ibeta_invb_f(float a, float x, float p)       { return boost::math::ibeta_invb(a, x, p); }

long double bs_beta_l(long double a, long double b)                      { return boost::math::beta(a, b); }
long double bs_fullBeta_l(long double a, long double b, long double x)   { return boost::math::beta(a, b, x); }
long double bs_ibeta_l(long double a, long double b, long double x)      { return boost::math::ibeta(a, b, x); }
long double bs_ibetac_l(long double a, long double b, long double x)     { return boost::math::ibetac(a, b, x); }
long double bs_ibeta_inv_l(long double a, long double b, long double p)  { return boost::math::ibeta_inv(a, b, p); }
long double bs_ibetac_inv_l(long double a, long double b, long double p) { return boost::math::ibetac_inv(a, b, p); }
long double bs_ibeta_inva_l(long double b, long double x, long double p) { return boost::math::ibeta_inva(b, x, p); }
long double bs_ibeta_invb_l(long double a, long double x, long double p) { return boost::math::ibeta_invb(a, x, p); }

// Digamma / Polygamma / Zeta
double bs_digamma(double x)            { return boost::math::digamma(x); }
double bs_trigamma(double x)           { return boost::math::trigamma(x); }
double bs_polygamma(int n, double x)   { return boost::math::polygamma(n, x); }
double bs_riemann_zeta(double x)       { return boost::math::zeta(x); }

float bs_digamma_f(float x)            { return boost::math::digamma(x); }
float bs_trigamma_f(float x)           { return boost::math::trigamma(x); }
float bs_polygamma_f(int n, float x)   { return boost::math::polygamma(n, x); }
float bs_riemann_zeta_f(float x)       { return boost::math::zeta(x); }

long double bs_digamma_l(long double x)          { return boost::math::digamma(x); }
long double bs_trigamma_l(long double x)         { return boost::math::trigamma(x); }
long double bs_polygamma_l(int n, long double x) { return boost::math::polygamma(n, x); }
long double bs_riemann_zeta_l(long double x)     { return boost::math::zeta(x); }

// Owen's T
double bs_owens_t(double h, double a)           { return boost::math::owens_t(h, a); }
float bs_owens_t_f(float h, float a)            { return boost::math::owens_t(h, a); }
long double bs_owens_t_l(long double h, long double a) { return boost::math::owens_t(h, a); }

// Exponential integrals and related
double bs_expint_En(int n, double x)   { return boost::math::expint(n, x); }
double bs_expm1(double x)              { return boost::math::expm1(x); }
double bs_log1p(double x)              { return boost::math::log1p(x); }
double bs_log1pmx(double x)            { return boost::math::log1pmx(x); }
double bs_powm1(double x, double y)    { return boost::math::powm1(x, y); }
double bs_cbrt(double x)               { return boost::math::cbrt(x); }

float bs_expint_En_f(int n, float x)   { return boost::math::expint(n, x); }
float bs_expm1_f(float x)              { return boost::math::expm1(x); }
float bs_log1p_f(float x)              { return boost::math::log1p(x); }
float bs_log1pmx_f(float x)            { return boost::math::log1pmx(x); }
float bs_powm1_f(float x, float y)     { return boost::math::powm1(x, y); }
float bs_cbrt_f(float x)               { return boost::math::cbrt(x); }

long double bs_expint_En_l(int n, long double x)   { return boost::math::expint(n, x); }
long double bs_expm1_l(long double x)              { return boost::math::expm1(x); }
long double bs_log1p_l(long double x)              { return boost::math::log1p(x); }
long double bs_log1pmx_l(long double x)            { return boost::math::log1pmx(x); }
long double bs_powm1_l(long double x, long double y) { return boost::math::powm1(x, y); }
long double bs_cbrt_l(long double x)               { return boost::math::cbrt(x); }

// Trig helpers
double bs_sin_pi(double x)             { return boost::math::sin_pi(x); }
double bs_cos_pi(double x)             { return boost::math::cos_pi(x); }

float bs_sin_pi_f(float x)             { return boost::math::sin_pi(x); }
float bs_cos_pi_f(float x)             { return boost::math::cos_pi(x); }

long double bs_sin_pi_l(long double x) { return boost::math::sin_pi(x); }
long double bs_cos_pi_l(long double x) { return boost::math::cos_pi(x); }

// Airy
double bs_airy_ai(double x)            { return boost::math::airy_ai(x); }
double bs_airy_bi(double x)            { return boost::math::airy_bi(x); }
double bs_airy_ai_prime(double x)      { return boost::math::airy_ai_prime(x); }
double bs_airy_bi_prime(double x)      { return boost::math::airy_bi_prime(x); }

float bs_airy_ai_f(float x)            { return boost::math::airy_ai(x); }
float bs_airy_bi_f(float x)            { return boost::math::airy_bi(x); }
float bs_airy_ai_prime_f(float x)      { return boost::math::airy_ai_prime(x); }
float bs_airy_bi_prime_f(float x)      { return boost::math::airy_bi_prime(x); }

long double bs_airy_ai_l(long double x)       { return boost::math::airy_ai(x); }
long double bs_airy_bi_l(long double x)       { return boost::math::airy_bi(x); }
long double bs_airy_ai_prime_l(long double x) { return boost::math::airy_ai_prime(x); }
long double bs_airy_bi_prime_l(long double x) { return boost::math::airy_bi_prime(x); }

// Bessel (cylindrical, real)
double bs_cyl_bessel_j(double v, double x) { return boost::math::cyl_bessel_j(v, x); }
double bs_cyl_neumann(double v, double x)  { return boost::math::cyl_neumann(v, x); }
double bs_cyl_bessel_i(double v, double x) { return boost::math::cyl_bessel_i(v, x); }
double bs_cyl_bessel_k(double v, double x) { return boost::math::cyl_bessel_k(v, x); }

float bs_cyl_bessel_j_f(float v, float x) { return boost::math::cyl_bessel_j(v, x); }
float bs_cyl_neumann_f(float v, float x)  { return boost::math::cyl_neumann(v, x); }
float bs_cyl_bessel_i_f(float v, float x) { return boost::math::cyl_bessel_i(v, x); }
float bs_cyl_bessel_k_f(float v, float x) { return boost::math::cyl_bessel_k(v, x); }

long double bs_cyl_bessel_j_l(long double v, long double x) { return boost::math::cyl_bessel_j(v, x); }
long double bs_cyl_neumann_l(long double v, long double x)  { return boost::math::cyl_neumann(v, x); }
long double bs_cyl_bessel_i_l(long double v, long double x) { return boost::math::cyl_bessel_i(v, x); }
long double bs_cyl_bessel_k_l(long double v, long double x) { return boost::math::cyl_bessel_k(v, x); }

// Legendre
double bs_legendre_p(int n, double x)                 { return boost::math::legendre_p(n, x); }
double bs_assoc_legendre_p(int n, int m, double x)    { return boost::math::legendre_p(n, m, x); }

float bs_legendre_p_f(int n, float x)                 { return boost::math::legendre_p(n, x); }
float bs_assoc_legendre_p_f(int n, int m, float x)    { return boost::math::legendre_p(n, m, x); }

long double bs_legendre_p_l(int n, long double x)              { return boost::math::legendre_p(n, x); }
long double bs_assoc_legendre_p_l(int n, int m, long double x) { return boost::math::legendre_p(n, m, x); }

// Elliptic integrals (Legendre forms)
double bs_ellint_1_complete(double k)                 { return boost::math::ellint_1(k); }
double bs_ellint_1(double k, double phi)              { return boost::math::ellint_1(k, phi); }
double bs_ellint_2_complete(double k)                 { return boost::math::ellint_2(k); }
double bs_ellint_2(double k, double phi)              { return boost::math::ellint_2(k, phi); }
double bs_ellint_3_complete(double k, double nu)      { return boost::math::ellint_3(k, nu); }
double bs_ellint_3(double k, double nu, double phi)   { return boost::math::ellint_3(k, nu, phi); }

float bs_ellint_1_complete_f(float k)                 { return boost::math::ellint_1(k); }
float bs_ellint_1_f(float k, float phi)               { return boost::math::ellint_1(k, phi); }
float bs_ellint_2_complete_f(float k)                 { return boost::math::ellint_2(k); }
float bs_ellint_2_f(float k, float phi)               { return boost::math::ellint_2(k, phi); }
float bs_ellint_3_complete_f(float k, float nu)       { return boost::math::ellint_3(k, nu); }
float bs_ellint_3_f(float k, float nu, float phi)     { return boost::math::ellint_3(k, nu, phi); }

long double bs_ellint_1_complete_l(long double k)                 { return boost::math::ellint_1(k); }
long double bs_ellint_1_l(long double k, long double phi)         { return boost::math::ellint_1(k, phi); }
long double bs_ellint_2_complete_l(long double k)                 { return boost::math::ellint_2(k); }
long double bs_ellint_2_l(long double k, long double phi)         { return boost::math::ellint_2(k, phi); }
long double bs_ellint_3_complete_l(long double k, long double nu) { return boost::math::ellint_3(k, nu); }
long double bs_ellint_3_l(long double k, long double nu, long double phi) { return boost::math::ellint_3(k, nu, phi); }

// Elliptic integrals (Carlson symmetric forms)
double bs_ellint_rc(double x, double y)               { return boost::math::ellint_rc(x, y); }
double bs_ellint_rf(double x, double y, double z)     { return boost::math::ellint_rf(x, y, z); }
double bs_ellint_rd(double x, double y, double z)     { return boost::math::ellint_rd(x, y, z); }
double bs_ellint_rj(double x, double y, double z, double p) { return boost::math::ellint_rj(x, y, z, p); }
double bs_ellint_rg(double x, double y, double z)     { return boost::math::ellint_rg(x, y, z); }

float bs_ellint_rc_f(float x, float y)                { return boost::math::ellint_rc(x, y); }
float bs_ellint_rf_f(float x, float y, float z)       { return boost::math::ellint_rf(x, y, z); }
float bs_ellint_rd_f(float x, float y, float z)       { return boost::math::ellint_rd(x, y, z); }
float bs_ellint_rj_f(float x, float y, float z, float p) { return boost::math::ellint_rj(x, y, z, p); }
float bs_ellint_rg_f(float x, float y, float z)       { return boost::math::ellint_rg(x, y, z); }

long double bs_ellint_rc_l(long double x, long double y)                { return boost::math::ellint_rc(x, y); }
long double bs_ellint_rf_l(long double x, long double y, long double z) { return boost::math::ellint_rf(x, y, z); }
long double bs_ellint_rd_l(long double x, long double y, long double z) { return boost::math::ellint_rd(x, y, z); }
long double bs_ellint_rj_l(long double x, long double y, long double z, long double p) { return boost::math::ellint_rj(x, y, z, p); }
long double bs_ellint_rg_l(long double x, long double y, long double z) { return boost::math::ellint_rg(x, y, z); }

// Lambert W (real branches)
double bs_lambert_w0(double x)        { return boost::math::lambert_w0(x); }
double bs_lambert_wm1(double x)       { return boost::math::lambert_wm1(x); }

float bs_lambert_w0_f(float x)        { return boost::math::lambert_w0(x); }
float bs_lambert_wm1_f(float x)       { return boost::math::lambert_wm1(x); }

long double bs_lambert_w0_l(long double x)  { return boost::math::lambert_w0(x); }
long double bs_lambert_wm1_l(long double x) { return boost::math::lambert_wm1(x); }

} // extern "C"
