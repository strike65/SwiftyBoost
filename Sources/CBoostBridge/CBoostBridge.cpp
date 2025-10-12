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
#include <boost/math/special_functions/hypergeometric_1F0.hpp>
#include <boost/math/special_functions/hypergeometric_0F1.hpp>
#include <boost/math/special_functions/hypergeometric_2F0.hpp>
#include <boost/math/special_functions/hypergeometric_1F1.hpp>
#include <boost/math/special_functions/hypergeometric_pFq.hpp>
#include <boost/math/special_functions/bernoulli.hpp>
#include <boost/math/special_functions/fibonacci.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/prime.hpp>

#include <vector>
#include <limits>
#include <stdexcept>

// Exception-safe wrapper: translate Boost throws to numeric sentinels.
// - Overflow -> +infinity
// - Domain/other errors -> quiet NaN
template <class T, class F>
inline T bs_wrap(F&& f) noexcept {
    try {
        return static_cast<T>(f());
    } catch (const std::overflow_error&) {
        return std::numeric_limits<T>::infinity();
    } catch (const std::domain_error&) {
        return std::numeric_limits<T>::quiet_NaN();
    } catch (...) {
        return std::numeric_limits<T>::quiet_NaN();
    }
}

extern "C" {

// Boost constants (double)
double bs_const_e(void)                 { return bs_wrap<double>([] { return boost::math::constants::e<double>(); }); }
double bs_const_pi(void)                { return bs_wrap<double>([] { return boost::math::constants::pi<double>(); }); }
double bs_const_two_pi(void)            { return bs_wrap<double>([] { return boost::math::constants::two_pi<double>(); }); }
double bs_const_half_pi(void)           { return bs_wrap<double>([] { return boost::math::constants::half_pi<double>(); }); }
double bs_const_third_pi(void)          { return bs_wrap<double>([] { return boost::math::constants::third_pi<double>(); }); }
double bs_const_two_thirds_pi(void)     { return bs_wrap<double>([] { return boost::math::constants::two_thirds_pi<double>(); }); }
double bs_const_quarter_pi(void)        { return bs_wrap<double>([] { return boost::math::constants::quarter_pi<double>(); }); }
double bs_const_three_quarters_pi(void) { return bs_wrap<double>([] { return boost::math::constants::three_quarters_pi<double>(); }); }
double bs_const_sixth_pi(void)          { return bs_wrap<double>([] { return boost::math::constants::sixth_pi<double>(); }); }
double bs_const_pi_sqr(void)            { return bs_wrap<double>([] { return boost::math::constants::pi_sqr<double>(); }); }

double bs_const_root_two(void)          { return bs_wrap<double>([] { return boost::math::constants::root_two<double>(); }); }
double bs_const_root_three(void)        { return bs_wrap<double>([] { return boost::math::constants::root_three<double>(); }); }
double bs_const_root_pi(void)           { return bs_wrap<double>([] { return boost::math::constants::root_pi<double>(); }); }

double bs_const_one_div_pi(void)        { return bs_wrap<double>([] { return boost::math::constants::one_div_pi<double>(); }); }
double bs_const_one_div_two_pi(void)    { return bs_wrap<double>([] { return boost::math::constants::one_div_two_pi<double>(); }); }
double bs_const_one_div_root_pi(void)   { return bs_wrap<double>([] { return boost::math::constants::one_div_root_pi<double>(); }); }
double bs_const_two_div_pi(void)        { return bs_wrap<double>([] { return boost::math::constants::two_div_pi<double>(); }); }
double bs_const_two_div_root_pi(void)   { return bs_wrap<double>([] { return boost::math::constants::two_div_root_pi<double>(); }); }

double bs_const_ln_two(void)            { return bs_wrap<double>([] { return boost::math::constants::ln_two<double>(); }); }
double bs_const_ln_ten(void)            { return bs_wrap<double>([] { return boost::math::constants::ln_ten<double>(); }); }
double bs_const_ln_ln_two(void)         { return bs_wrap<double>([] { return boost::math::constants::ln_ln_two<double>(); }); }

double bs_const_euler(void)             { return bs_wrap<double>([] { return boost::math::constants::euler<double>(); }); }
double bs_const_catalan(void)           { return bs_wrap<double>([] { return boost::math::constants::catalan<double>(); }); }
double bs_const_zeta_three(void)        { return bs_wrap<double>([] { return boost::math::constants::zeta_three<double>(); }); }
double bs_const_phi(void)               { return bs_wrap<double>([] { return boost::math::constants::phi<double>(); }); }

// Boost constants (float)
float bs_const_e_f(void)                 { return bs_wrap<float>([] { return boost::math::constants::e<float>(); }); }
float bs_const_pi_f(void)                { return bs_wrap<float>([] { return boost::math::constants::pi<float>(); }); }
float bs_const_two_pi_f(void)            { return bs_wrap<float>([] { return boost::math::constants::two_pi<float>(); }); }
float bs_const_half_pi_f(void)           { return bs_wrap<float>([] { return boost::math::constants::half_pi<float>(); }); }
float bs_const_third_pi_f(void)          { return bs_wrap<float>([] { return boost::math::constants::third_pi<float>(); }); }
float bs_const_two_thirds_pi_f(void)     { return bs_wrap<float>([] { return boost::math::constants::two_thirds_pi<float>(); }); }
float bs_const_quarter_pi_f(void)        { return bs_wrap<float>([] { return boost::math::constants::quarter_pi<float>(); }); }
float bs_const_three_quarters_pi_f(void) { return bs_wrap<float>([] { return boost::math::constants::three_quarters_pi<float>(); }); }
float bs_const_sixth_pi_f(void)          { return bs_wrap<float>([] { return boost::math::constants::sixth_pi<float>(); }); }
float bs_const_pi_sqr_f(void)            { return bs_wrap<float>([] { return boost::math::constants::pi_sqr<float>(); }); }

float bs_const_root_two_f(void)          { return bs_wrap<float>([] { return boost::math::constants::root_two<float>(); }); }
float bs_const_root_three_f(void)        { return bs_wrap<float>([] { return boost::math::constants::root_three<float>(); }); }
float bs_const_root_pi_f(void)           { return bs_wrap<float>([] { return boost::math::constants::root_pi<float>(); }); }

float bs_const_one_div_pi_f(void)        { return bs_wrap<float>([] { return boost::math::constants::one_div_pi<float>(); }); }
float bs_const_one_div_two_pi_f(void)    { return bs_wrap<float>([] { return boost::math::constants::one_div_two_pi<float>(); }); }
float bs_const_one_div_root_pi_f(void)   { return bs_wrap<float>([] { return boost::math::constants::one_div_root_pi<float>(); }); }
float bs_const_two_div_pi_f(void)        { return bs_wrap<float>([] { return boost::math::constants::two_div_pi<float>(); }); }
float bs_const_two_div_root_pi_f(void)   { return bs_wrap<float>([] { return boost::math::constants::two_div_root_pi<float>(); }); }

float bs_const_ln_two_f(void)            { return bs_wrap<float>([] { return boost::math::constants::ln_two<float>(); }); }
float bs_const_ln_ten_f(void)            { return bs_wrap<float>([] { return boost::math::constants::ln_ten<float>(); }); }
float bs_const_ln_ln_two_f(void)         { return bs_wrap<float>([] { return boost::math::constants::ln_ln_two<float>(); }); }

float bs_const_euler_f(void)             { return bs_wrap<float>([] { return boost::math::constants::euler<float>(); }); }
float bs_const_catalan_f(void)           { return bs_wrap<float>([] { return boost::math::constants::catalan<float>(); }); }
float bs_const_zeta_three_f(void)        { return bs_wrap<float>([] { return boost::math::constants::zeta_three<float>(); }); }
float bs_const_phi_f(void)               { return bs_wrap<float>([] { return boost::math::constants::phi<float>(); }); }

// Boost constants (long double)
long double bs_const_e_l(void)                 { return bs_wrap<long double>([] { return boost::math::constants::e<long double>(); }); }
long double bs_const_pi_l(void)                { return bs_wrap<long double>([] { return boost::math::constants::pi<long double>(); }); }
long double bs_const_two_pi_l(void)            { return bs_wrap<long double>([] { return boost::math::constants::two_pi<long double>(); }); }
long double bs_const_half_pi_l(void)           { return bs_wrap<long double>([] { return boost::math::constants::half_pi<long double>(); }); }
long double bs_const_third_pi_l(void)          { return bs_wrap<long double>([] { return boost::math::constants::third_pi<long double>(); }); }
long double bs_const_two_thirds_pi_l(void)     { return bs_wrap<long double>([] { return boost::math::constants::two_thirds_pi<long double>(); }); }
long double bs_const_quarter_pi_l(void)        { return bs_wrap<long double>([] { return boost::math::constants::quarter_pi<long double>(); }); }
long double bs_const_three_quarters_pi_l(void) { return bs_wrap<long double>([] { return boost::math::constants::three_quarters_pi<long double>(); }); }
long double bs_const_sixth_pi_l(void)          { return bs_wrap<long double>([] { return boost::math::constants::sixth_pi<long double>(); }); }
long double bs_const_pi_sqr_l(void)            { return bs_wrap<long double>([] { return boost::math::constants::pi_sqr<long double>(); }); }

long double bs_const_root_two_l(void)          { return bs_wrap<long double>([] { return boost::math::constants::root_two<long double>(); }); }
long double bs_const_root_three_l(void)        { return bs_wrap<long double>([] { return boost::math::constants::root_three<long double>(); }); }
long double bs_const_root_pi_l(void)           { return bs_wrap<long double>([] { return boost::math::constants::root_pi<long double>(); }); }

long double bs_const_one_div_pi_l(void)        { return bs_wrap<long double>([] { return boost::math::constants::one_div_pi<long double>(); }); }
long double bs_const_one_div_two_pi_l(void)    { return bs_wrap<long double>([] { return boost::math::constants::one_div_two_pi<long double>(); }); }
long double bs_const_one_div_root_pi_l(void)   { return bs_wrap<long double>([] { return boost::math::constants::one_div_root_pi<long double>(); }); }
long double bs_const_two_div_pi_l(void)        { return bs_wrap<long double>([] { return boost::math::constants::two_div_pi<long double>(); }); }
long double bs_const_two_div_root_pi_l(void)   { return bs_wrap<long double>([] { return boost::math::constants::two_div_root_pi<long double>(); }); }

long double bs_const_ln_two_l(void)            { return bs_wrap<long double>([] { return boost::math::constants::ln_two<long double>(); }); }
long double bs_const_ln_ten_l(void)            { return bs_wrap<long double>([] { return boost::math::constants::ln_ten<long double>(); }); }
long double bs_const_ln_ln_two_l(void)         { return bs_wrap<long double>([] { return boost::math::constants::ln_ln_two<long double>(); }); }

long double bs_const_euler_l(void)             { return bs_wrap<long double>([] { return boost::math::constants::euler<long double>(); }); }
long double bs_const_catalan_l(void)           { return bs_wrap<long double>([] { return boost::math::constants::catalan<long double>(); }); }
long double bs_const_zeta_three_l(void)        { return bs_wrap<long double>([] { return boost::math::constants::zeta_three<long double>(); }); }
long double bs_const_phi_l(void)               { return bs_wrap<long double>([] { return boost::math::constants::phi<long double>(); }); }

// Gamma / Error
double bs_tgamma(double x)              { return bs_wrap<double>([&] { return boost::math::tgamma(x); }); }
double bs_lgamma(double x)              { return bs_wrap<double>([&] { return boost::math::lgamma(x); }); }
double bs_erf(double x)                 { return bs_wrap<double>([&] { return boost::math::erf(x); }); }
double bs_erfc(double x)                { return bs_wrap<double>([&] { return boost::math::erfc(x); }); }

float bs_tgamma_f(float x)              { return bs_wrap<float>([&] { return boost::math::tgamma(x); }); }
float bs_lgamma_f(float x)              { return bs_wrap<float>([&] { return boost::math::lgamma(x); }); }
float bs_erf_f(float x)                 { return bs_wrap<float>([&] { return boost::math::erf(x); }); }
float bs_erfc_f(float x)                { return bs_wrap<float>([&] { return boost::math::erfc(x); }); }

long double bs_tgamma_l(long double x)  { return bs_wrap<long double>([&] { return boost::math::tgamma(x); }); }
long double bs_lgamma_l(long double x)  { return bs_wrap<long double>([&] { return boost::math::lgamma(x); }); }
long double bs_erf_l(long double x)     { return bs_wrap<long double>([&] { return boost::math::erf(x); }); }
long double bs_erfc_l(long double x)    { return bs_wrap<long double>([&] { return boost::math::erfc(x); }); }

// Incomplete gamma (lower/upper, regularized, and inverses)
double bs_tgamma_lower(double a, double x)      { return bs_wrap<double>([&] { return boost::math::tgamma_lower(a, x); }); }
double bs_tgamma_upper(double a, double x)      { return bs_wrap<double>([&] { return boost::math::tgamma(a, x); }); }
double bs_gamma_p(double a, double x)           { return bs_wrap<double>([&] { return boost::math::gamma_p(a, x); }); }
double bs_gamma_q(double a, double x)           { return bs_wrap<double>([&] { return boost::math::gamma_q(a, x); }); }
double bs_gamma_p_inv(double a, double p)       { return bs_wrap<double>([&] { return boost::math::gamma_p_inv(a, p); }); }
double bs_gamma_q_inv(double a, double q)       { return bs_wrap<double>([&] { return boost::math::gamma_q_inv(a, q); }); }

// Derivatives of regularized incomplete gamma functions
double bs_gamma_p_derivative(double a, double x) { return bs_wrap<double>([&] { return boost::math::gamma_p_derivative(a, x); }); }

float bs_tgamma_lower_f(float a, float x)       { return bs_wrap<float>([&] { return boost::math::tgamma_lower(a, x); }); }
float bs_tgamma_upper_f(float a, float x)       { return bs_wrap<float>([&] { return boost::math::tgamma(a, x); }); }
float bs_gamma_p_f(float a, float x)            { return bs_wrap<float>([&] { return boost::math::gamma_p(a, x); }); }
float bs_gamma_q_f(float a, float x)            { return bs_wrap<float>([&] { return boost::math::gamma_q(a, x); }); }
float bs_gamma_p_inv_f(float a, float p)        { return bs_wrap<float>([&] { return boost::math::gamma_p_inv(a, p); }); }
float bs_gamma_q_inv_f(float a, float q)        { return bs_wrap<float>([&] { return boost::math::gamma_q_inv(a, q); }); }

// Derivatives (float)
float bs_gamma_p_derivative_f(float a, float x) { return bs_wrap<float>([&] { return boost::math::gamma_p_derivative(a, x); }); }

long double bs_tgamma_lower_l(long double a, long double x) { return bs_wrap<long double>([&] { return boost::math::tgamma_lower(a, x); }); }
long double bs_tgamma_upper_l(long double a, long double x) { return bs_wrap<long double>([&] { return boost::math::tgamma(a, x); }); }
long double bs_gamma_p_l(long double a, long double x)      { return bs_wrap<long double>([&] { return boost::math::gamma_p(a, x); }); }
long double bs_gamma_q_l(long double a, long double x)      { return bs_wrap<long double>([&] { return boost::math::gamma_q(a, x); }); }
long double bs_gamma_p_inv_l(long double a, long double p)  { return bs_wrap<long double>([&] { return boost::math::gamma_p_inv(a, p); }); }
long double bs_gamma_q_inv_l(long double a, long double q)  { return bs_wrap<long double>([&] { return boost::math::gamma_q_inv(a, q); }); }

// Derivatives (long double)
long double bs_gamma_p_derivative_l(long double a, long double x) { return bs_wrap<long double>([&] { return boost::math::gamma_p_derivative(a, x); }); }

// Beta family
double bs_beta(double a, double b)                     { return bs_wrap<double>([&] { return boost::math::beta(a, b); }); }
double bs_fullBeta(double a, double b, double x)       { return bs_wrap<double>([&] { return boost::math::beta(a, b, x); }); }
double bs_ibeta(double a, double b, double x)          { return bs_wrap<double>([&] { return boost::math::ibeta(a, b, x); }); }
double bs_ibetac(double a, double b, double x)         { return bs_wrap<double>([&] { return boost::math::ibetac(a, b, x); }); }
double bs_ibeta_inv(double a, double b, double p)      { return bs_wrap<double>([&] { return boost::math::ibeta_inv(a, b, p); }); }
double bs_ibetac_inv(double a, double b, double p)     { return bs_wrap<double>([&] { return boost::math::ibetac_inv(a, b, p); }); }
double bs_ibeta_inva(double b, double x, double p)     { return bs_wrap<double>([&] { return boost::math::ibeta_inva(b, x, p); }); }
double bs_ibeta_invb(double a, double x, double p)     { return bs_wrap<double>([&] { return boost::math::ibeta_invb(a, x, p); }); }
double bs_ibeta_derivative(double a, double b, double x) { return bs_wrap<double>([&] { return boost::math::ibeta_derivative(a, b, x); }); }

float bs_beta_f(float a, float b)                      { return bs_wrap<float>([&] { return boost::math::beta(a, b); }); }
float bs_fullBeta_f(float a, float b, float x)         { return bs_wrap<float>([&] { return boost::math::beta(a, b, x); }); }
float bs_ibeta_f(float a, float b, float x)            { return bs_wrap<float>([&] { return boost::math::ibeta(a, b, x); }); }
float bs_ibetac_f(float a, float b, float x)           { return bs_wrap<float>([&] { return boost::math::ibetac(a, b, x); }); }
float bs_ibeta_inv_f(float a, float b, float p)        { return bs_wrap<float>([&] { return boost::math::ibeta_inv(a, b, p); }); }
float bs_ibetac_inv_f(float a, float b, float p)       { return bs_wrap<float>([&] { return boost::math::ibetac_inv(a, b, p); }); }
float bs_ibeta_inva_f(float b, float x, float p)       { return bs_wrap<float>([&] { return boost::math::ibeta_inva(b, x, p); }); }
float bs_ibeta_invb_f(float a, float x, float p)       { return bs_wrap<float>([&] { return boost::math::ibeta_invb(a, x, p); }); }
float bs_ibeta_derivative_f(float a, float b, float x) { return bs_wrap<float>([&] { return boost::math::ibeta_derivative(a, b, x); }); }

long double bs_beta_l(long double a, long double b)                      { return bs_wrap<long double>([&] { return boost::math::beta(a, b); }); }
long double bs_fullBeta_l(long double a, long double b, long double x)   { return bs_wrap<long double>([&] { return boost::math::beta(a, b, x); }); }
long double bs_ibeta_l(long double a, long double b, long double x)      { return bs_wrap<long double>([&] { return boost::math::ibeta(a, b, x); }); }
long double bs_ibetac_l(long double a, long double b, long double x)     { return bs_wrap<long double>([&] { return boost::math::ibetac(a, b, x); }); }
long double bs_ibeta_inv_l(long double a, long double b, long double p)  { return bs_wrap<long double>([&] { return boost::math::ibeta_inv(a, b, p); }); }
long double bs_ibetac_inv_l(long double a, long double b, long double p) { return bs_wrap<long double>([&] { return boost::math::ibetac_inv(a, b, p); }); }
long double bs_ibeta_inva_l(long double b, long double x, long double p) { return bs_wrap<long double>([&] { return boost::math::ibeta_inva(b, x, p); }); }
long double bs_ibeta_invb_l(long double a, long double x, long double p) { return bs_wrap<long double>([&] { return boost::math::ibeta_invb(a, x, p); }); }
long double bs_ibeta_derivative_l(long double a, long double b, long double x) { return bs_wrap<long double>([&] { return boost::math::ibeta_derivative(a, b, x); }); }

// Digamma / Polygamma / Zeta
double bs_digamma(double x)            { return bs_wrap<double>([&] { return boost::math::digamma(x); }); }
double bs_trigamma(double x)           { return bs_wrap<double>([&] { return boost::math::trigamma(x); }); }
double bs_polygamma(int n, double x)   { return bs_wrap<double>([&] { return boost::math::polygamma(n, x); }); }
double bs_riemann_zeta(double x)       { return bs_wrap<double>([&] { return boost::math::zeta(x); }); }

float bs_digamma_f(float x)            { return bs_wrap<float>([&] { return boost::math::digamma(x); }); }
float bs_trigamma_f(float x)           { return bs_wrap<float>([&] { return boost::math::trigamma(x); }); }
float bs_polygamma_f(int n, float x)   { return bs_wrap<float>([&] { return boost::math::polygamma(n, x); }); }
float bs_riemann_zeta_f(float x)       { return bs_wrap<float>([&] { return boost::math::zeta(x); }); }

long double bs_digamma_l(long double x)          { return bs_wrap<long double>([&] { return boost::math::digamma(x); }); }
long double bs_trigamma_l(long double x)         { return bs_wrap<long double>([&] { return boost::math::trigamma(x); }); }
long double bs_polygamma_l(int n, long double x) { return bs_wrap<long double>([&] { return boost::math::polygamma(n, x); }); }
long double bs_riemann_zeta_l(long double x)     { return bs_wrap<long double>([&] { return boost::math::zeta(x); }); }

// Owen's T
double bs_owens_t(double h, double a)           { return bs_wrap<double>([&] { return boost::math::owens_t(h, a); }); }
float bs_owens_t_f(float h, float a)            { return bs_wrap<float>([&] { return boost::math::owens_t(h, a); }); }
long double bs_owens_t_l(long double h, long double a) { return bs_wrap<long double>([&] { return boost::math::owens_t(h, a); }); }

// Exponential integrals and related
double bs_expint_Ei(double x)        { return bs_wrap<double>([&] { return boost::math::expint(x); }); }
double bs_expint_En(int n, double x) { return bs_wrap<double>([&] { return boost::math::expint(n, x); }); }
double bs_expm1(double x)            { return bs_wrap<double>([&] { return boost::math::expm1(x); }); }
double bs_log1p(double x)            { return bs_wrap<double>([&] { return boost::math::log1p(x); }); }
double bs_log1pmx(double x)          { return bs_wrap<double>([&] { return boost::math::log1pmx(x); }); }
double bs_powm1(double x, double y)  { return bs_wrap<double>([&] { return boost::math::powm1(x, y); }); }
double bs_cbrt(double x)             { return bs_wrap<double>([&] { return boost::math::cbrt(x); }); }

float bs_expint_Ei_f(float x)        { return bs_wrap<float>([&] { return boost::math::expint(x); }); }
float bs_expint_En_f(int n, float x) { return bs_wrap<float>([&] { return boost::math::expint(n, x); }); }
float bs_expm1_f(float x)            { return bs_wrap<float>([&] { return boost::math::expm1(x); }); }
float bs_log1p_f(float x)            { return bs_wrap<float>([&] { return boost::math::log1p(x); }); }
float bs_log1pmx_f(float x)          { return bs_wrap<float>([&] { return boost::math::log1pmx(x); }); }
float bs_powm1_f(float x, float y)   { return bs_wrap<float>([&] { return boost::math::powm1(x, y); }); }
float bs_cbrt_f(float x)             { return bs_wrap<float>([&] { return boost::math::cbrt(x); }); }

long double bs_expint_Ei_l(long double x)          { return bs_wrap<long double>([&] { return boost::math::expint(x); }); }
long double bs_expint_En_l(int n, long double x)   { return bs_wrap<long double>([&] { return boost::math::expint(n, x); }); }
long double bs_expm1_l(long double x)              { return bs_wrap<long double>([&] { return boost::math::expm1(x); }); }
long double bs_log1p_l(long double x)              { return bs_wrap<long double>([&] { return boost::math::log1p(x); }); }
long double bs_log1pmx_l(long double x)            { return bs_wrap<long double>([&] { return boost::math::log1pmx(x); }); }
long double bs_powm1_l(long double x, long double y) { return bs_wrap<long double>([&] { return boost::math::powm1(x, y); }); }
long double bs_cbrt_l(long double x)               { return bs_wrap<long double>([&] { return boost::math::cbrt(x); }); }

// Trig helpers
double bs_sin_pi(double x)             { return bs_wrap<double>([&] { return boost::math::sin_pi(x); }); }
double bs_cos_pi(double x)             { return bs_wrap<double>([&] { return boost::math::cos_pi(x); }); }

float bs_sin_pi_f(float x)             { return bs_wrap<float>([&] { return boost::math::sin_pi(x); }); }
float bs_cos_pi_f(float x)             { return bs_wrap<float>([&] { return boost::math::cos_pi(x); }); }

long double bs_sin_pi_l(long double x) { return bs_wrap<long double>([&] { return boost::math::sin_pi(x); }); }
long double bs_cos_pi_l(long double x) { return bs_wrap<long double>([&] { return boost::math::cos_pi(x); }); }

// Airy
double bs_airy_ai(double x)            { return bs_wrap<double>([&] { return boost::math::airy_ai(x); }); }
double bs_airy_bi(double x)            { return bs_wrap<double>([&] { return boost::math::airy_bi(x); }); }
double bs_airy_ai_prime(double x)      { return bs_wrap<double>([&] { return boost::math::airy_ai_prime(x); }); }
double bs_airy_bi_prime(double x)      { return bs_wrap<double>([&] { return boost::math::airy_bi_prime(x); }); }

float bs_airy_ai_f(float x)            { return bs_wrap<float>([&] { return boost::math::airy_ai(x); }); }
float bs_airy_bi_f(float x)            { return bs_wrap<float>([&] { return boost::math::airy_bi(x); }); }
float bs_airy_ai_prime_f(float x)      { return bs_wrap<float>([&] { return boost::math::airy_ai_prime(x); }); }
float bs_airy_bi_prime_f(float x)      { return bs_wrap<float>([&] { return boost::math::airy_bi_prime(x); }); }

long double bs_airy_ai_l(long double x)       { return bs_wrap<long double>([&] { return boost::math::airy_ai(x); }); }
long double bs_airy_bi_l(long double x)       { return bs_wrap<long double>([&] { return boost::math::airy_bi(x); }); }
long double bs_airy_ai_prime_l(long double x) { return bs_wrap<long double>([&] { return boost::math::airy_ai_prime(x); }); }
long double bs_airy_bi_prime_l(long double x) { return bs_wrap<long double>([&] { return boost::math::airy_bi_prime(x); }); }

// Bessel (cylindrical, real)
double bs_cyl_bessel_j(double v, double x) { return bs_wrap<double>([&] { return boost::math::cyl_bessel_j(v, x); }); }
double bs_cyl_neumann(double v, double x)  { return bs_wrap<double>([&] { return boost::math::cyl_neumann(v, x); }); }
double bs_cyl_bessel_i(double v, double x) { return bs_wrap<double>([&] { return boost::math::cyl_bessel_i(v, x); }); }
double bs_cyl_bessel_k(double v, double x) { return bs_wrap<double>([&] { return boost::math::cyl_bessel_k(v, x); }); }

float bs_cyl_bessel_j_f(float v, float x) { return bs_wrap<float>([&] { return boost::math::cyl_bessel_j(v, x); }); }
float bs_cyl_neumann_f(float v, float x)  { return bs_wrap<float>([&] { return boost::math::cyl_neumann(v, x); }); }
float bs_cyl_bessel_i_f(float v, float x) { return bs_wrap<float>([&] { return boost::math::cyl_bessel_i(v, x); }); }
float bs_cyl_bessel_k_f(float v, float x) { return bs_wrap<float>([&] { return boost::math::cyl_bessel_k(v, x); }); }

long double bs_cyl_bessel_j_l(long double v, long double x) { return bs_wrap<long double>([&] { return boost::math::cyl_bessel_j(v, x); }); }
long double bs_cyl_neumann_l(long double v, long double x)  { return bs_wrap<long double>([&] { return boost::math::cyl_neumann(v, x); }); }
long double bs_cyl_bessel_i_l(long double v, long double x) { return bs_wrap<long double>([&] { return boost::math::cyl_bessel_i(v, x); }); }
long double bs_cyl_bessel_k_l(long double v, long double x) { return bs_wrap<long double>([&] { return boost::math::cyl_bessel_k(v, x); }); }

// Legendre
double bs_legendre_p(int n, double x)                 { return bs_wrap<double>([&] { return boost::math::legendre_p(n, x); }); }
double bs_assoc_legendre_p(int n, int m, double x)    { return bs_wrap<double>([&] { return boost::math::legendre_p(n, m, x); }); }

float bs_legendre_p_f(int n, float x)                 { return bs_wrap<float>([&] { return boost::math::legendre_p(n, x); }); }
float bs_assoc_legendre_p_f(int n, int m, float x)    { return bs_wrap<float>([&] { return boost::math::legendre_p(n, m, x); }); }

long double bs_legendre_p_l(int n, long double x)              { return bs_wrap<long double>([&] { return boost::math::legendre_p(n, x); }); }
long double bs_assoc_legendre_p_l(int n, int m, long double x) { return bs_wrap<long double>([&] { return boost::math::legendre_p(n, m, x); }); }

// Elliptic integrals (Legendre forms)
double bs_ellint_1_complete(double k)                 { return bs_wrap<double>([&] { return boost::math::ellint_1(k); }); }
double bs_ellint_1(double k, double phi)              { return bs_wrap<double>([&] { return boost::math::ellint_1(k, phi); }); }
double bs_ellint_2_complete(double k)                 { return bs_wrap<double>([&] { return boost::math::ellint_2(k); }); }
double bs_ellint_2(double k, double phi)              { return bs_wrap<double>([&] { return boost::math::ellint_2(k, phi); }); }
double bs_ellint_3(double k, double nu, double phi)   { return bs_wrap<double>([&] { return boost::math::ellint_3(k, nu, phi); }); } // incomplete elliptic integral of the third kind ∏(k, nu, phi)
double bs_ellint_3_complete(double k, double nu)      { return bs_wrap<double>([&] { return boost::math::ellint_3(k, nu); }); }


float bs_ellint_1_complete_f(float k)                 { return bs_wrap<float>([&] { return boost::math::ellint_1(k); }); }
float bs_ellint_1_f(float k, float phi)               { return bs_wrap<float>([&] { return boost::math::ellint_1(k, phi); }); }
float bs_ellint_2_complete_f(float k)                 { return bs_wrap<float>([&] { return boost::math::ellint_2(k); }); }
float bs_ellint_2_f(float k, float phi)               { return bs_wrap<float>([&] { return boost::math::ellint_2(k, phi); }); }
float bs_ellint_3_f(float k, float nu, float phi)     { return bs_wrap<float>([&] { return boost::math::ellint_3(k, nu, phi); }); } // incomplete elliptic integral of the third kind ∏(k, nu, phi)
float bs_ellint_3_complete_f(float k, float nu)      { return bs_wrap<float>([&] { return boost::math::ellint_3(k, nu); }); }

long double bs_ellint_1_complete_l(long double k)                 { return bs_wrap<long double>([&] { return boost::math::ellint_1(k); }); }
long double bs_ellint_1_l(long double k, long double phi)         { return bs_wrap<long double>([&] { return boost::math::ellint_1(k, phi); }); }
long double bs_ellint_2_complete_l(long double k)                 { return bs_wrap<long double>([&] { return boost::math::ellint_2(k); }); }
long double bs_ellint_2_l(long double k, long double phi)         { return bs_wrap<long double>([&] { return boost::math::ellint_2(k, phi); }); }
long double bs_ellint_3_l(long double k, long double nu,long double phi) { return bs_wrap<long double>([&] { return boost::math::ellint_3(k, nu, phi); }); } // incomplete elliptic integral of the third kind ∏(k, nu, phi)
long double bs_ellint_3_complete_l(long double k, long double nu)      { return bs_wrap<long double>([&] { return boost::math::ellint_3(k, nu); }); }

// Elliptic integrals (Carlson symmetric forms)
double bs_ellint_rc(double x, double y)               { return bs_wrap<double>([&] { return boost::math::ellint_rc(x, y); }); }
double bs_ellint_rf(double x, double y, double z)     { return bs_wrap<double>([&] { return boost::math::ellint_rf(x, y, z); }); }
double bs_ellint_rd(double x, double y, double z)     { return bs_wrap<double>([&] { return boost::math::ellint_rd(x, y, z); }); }
double bs_ellint_rj(double x, double y, double z, double p) { return bs_wrap<double>([&] { return boost::math::ellint_rj(x, y, z, p); }); }
double bs_ellint_rg(double x, double y, double z)     { return bs_wrap<double>([&] { return boost::math::ellint_rg(x, y, z); }); }

float bs_ellint_rc_f(float x, float y)                { return bs_wrap<float>([&] { return boost::math::ellint_rc(x, y); }); }
float bs_ellint_rf_f(float x, float y, float z)       { return bs_wrap<float>([&] { return boost::math::ellint_rf(x, y, z); }); }
float bs_ellint_rd_f(float x, float y, float z)       { return bs_wrap<float>([&] { return boost::math::ellint_rd(x, y, z); }); }
float bs_ellint_rj_f(float x, float y, float z, float p) { return bs_wrap<float>([&] { return boost::math::ellint_rj(x, y, z, p); }); }
float bs_ellint_rg_f(float x, float y, float z)       { return bs_wrap<float>([&] { return boost::math::ellint_rg(x, y, z); }); }

long double bs_ellint_rc_l(long double x, long double y)                { return bs_wrap<long double>([&] { return boost::math::ellint_rc(x, y); }); }
long double bs_ellint_rf_l(long double x, long double y, long double z) { return bs_wrap<long double>([&] { return boost::math::ellint_rf(x, y, z); }); }
long double bs_ellint_rd_l(long double x, long double y, long double z) { return bs_wrap<long double>([&] { return boost::math::ellint_rd(x, y, z); }); }
long double bs_ellint_rj_l(long double x, long double y, long double z, long double p) { return bs_wrap<long double>([&] { return boost::math::ellint_rj(x, y, z, p); }); }
long double bs_ellint_rg_l(long double x, long double y, long double z) { return bs_wrap<long double>([&] { return boost::math::ellint_rg(x, y, z); }); }

// Lambert W (real branches)
double bs_lambert_w0(double x)        { return bs_wrap<double>([&] { return boost::math::lambert_w0(x); }); }
double bs_lambert_wm1(double x)       { return bs_wrap<double>([&] { return boost::math::lambert_wm1(x); }); }

float bs_lambert_w0_f(float x)        { return bs_wrap<float>([&] { return boost::math::lambert_w0(x); }); }
float bs_lambert_wm1_f(float x)       { return bs_wrap<float>([&] { return boost::math::lambert_wm1(x); }); }

long double bs_lambert_w0_l(long double x)  { return bs_wrap<long double>([&] { return boost::math::lambert_w0(x); }); }
long double bs_lambert_wm1_l(long double x) { return bs_wrap<long double>([&] { return boost::math::lambert_wm1(x); }); }

// Γ(a) / Γ(b)
double bs_tgamma_ratio(double a, double b) { return bs_wrap<double>([&] { return boost::math::tgamma_ratio(a, b); }); }
float  bs_tgamma_ratio_f(float a, float b) { return bs_wrap<float>([&] { return boost::math::tgamma_ratio(a, b); }); }
long double bs_tgamma_ratio_l(long double a, long double b) { return bs_wrap<long double>([&] { return boost::math::tgamma_ratio(a, b); }); }

// Γ(a) / Γ(a + delta)
double bs_tgamma_delta_ratio(double a, double delta)  { return bs_wrap<double>([&] { return boost::math::tgamma_delta_ratio(a, delta); }); }
float  bs_tgamma_delta_ratio_f(float a, float delta)  { return bs_wrap<float>([&] { return boost::math::tgamma_delta_ratio(a, delta); }); }
long double bs_tgamma_delta_ratio_l(long double a, long double delta)  { return bs_wrap<long double>([&] { return boost::math::tgamma_delta_ratio(a, delta); }); }

// Hypergeometric functions
double bs_hypergeometric_1F0(double a, double z) { return bs_wrap<double>([&] { return boost::math::hypergeometric_1F0(a, z); }); }
float  bs_hypergeometric_1F0_f(float a, float z) { return bs_wrap<float>([&] { return boost::math::hypergeometric_1F0(a, z); }); }
long double bs_hypergeometric_1F0_l(long double a, long double z) { return bs_wrap<long double>([&] { return boost::math::hypergeometric_1F0(a, z); }); }

double bs_hypergeometric_0F1(double b, double z) { return bs_wrap<double>([&] { return boost::math::hypergeometric_0F1(b, z); }); }
float  bs_hypergeometric_0F1_f(float b, float z) { return bs_wrap<float>([&] { return boost::math::hypergeometric_0F1(b, z); }); }
long double bs_hypergeometric_0F1_l(long double b, long double z) { return bs_wrap<long double>([&] { return boost::math::hypergeometric_0F1(b, z); }); }

double bs_hypergeometric_2F0(double a, double b, double z) { return bs_wrap<double>([&] { return boost::math::hypergeometric_2F0(a, b, z); }); }
float  bs_hypergeometric_2F0_f(float a, float b, float z)  { return bs_wrap<float>([&] { return boost::math::hypergeometric_2F0(a, b, z); }); }
long double bs_hypergeometric_2F0_l(long double a, long double b, long double z) { return bs_wrap<long double>([&] { return boost::math::hypergeometric_2F0(a, b, z); }); }

double bs_hypergeometric_1F1(double a, double b, double z) { return bs_wrap<double>([&] { return boost::math::hypergeometric_1F1(a, b, z); }); }
float  bs_hypergeometric_1F1_f(float a, float b, float z)  { return bs_wrap<float>([&] { return boost::math::hypergeometric_1F1(a, b, z); }); }
long double bs_hypergeometric_1F1_l(long double a, long double b, long double z) { return bs_wrap<long double>([&] { return boost::math::hypergeometric_1F1(a, b, z); }); }

double bs_hypergeometric_pFq(const double* a, size_t p, const double* b, size_t q, double z) {
    return bs_wrap<double>([&] {
        std::vector<double> va(a, a + p);
        std::vector<double> vb(b, b + q);
        return boost::math::hypergeometric_pFq(va, vb, z);
    });
}
float bs_hypergeometric_pFq_f(const float* a, size_t p, const float* b, size_t q, float z) {
    return bs_wrap<float>([&] {
        std::vector<float> va(a, a + p);
        std::vector<float> vb(b, b + q);
        return boost::math::hypergeometric_pFq(va, vb, z);
    });
}
long double bs_hypergeometric_pFq_l(const long double* a, size_t p, const long double* b, size_t q, long double z) {
    return bs_wrap<long double>([&] {
        std::vector<long double> va(a, a + p);
        std::vector<long double> vb(b, b + q);
        return boost::math::hypergeometric_pFq(va, vb, z);
    });
}

// Bernoulli Numbers
double bs_bernoulli_b2n(const int n) {
    return bs_wrap<double>([&] { return boost::math::bernoulli_b2n<double>(n); });
}
float bs_bernoulli_b2n_f(const int n) {
    return bs_wrap<float>([&] { return boost::math::bernoulli_b2n<float>(n); });
}
long double bs_bernoulli_b2n_l(const int n) {
    return bs_wrap<long double>([&] { return boost::math::bernoulli_b2n<long double>(n); });
}

// Tangent Numbers (scalar)
double bs_tangent_t2n(const int n) {
    return bs_wrap<double>([&] { return boost::math::tangent_t2n<double>(n); });
}
float bs_tangent_t2n_f(const int n) {
    return bs_wrap<float>([&] { return boost::math::tangent_t2n<float>(n); });
}
long double bs_tangent_t2n_l(const int n) {
    return bs_wrap<long double>([&] { return boost::math::tangent_t2n<long double>(n); });
}

// Tangent Numbers (bulk sequence)
void bs_tangent_t2n_seq(int start_index, unsigned int count, double* out) {
    if (!out || count == 0) return;
    try {
        boost::math::tangent_t2n<double>(start_index, count, out);
    } catch (...) {
        // Leave output as-is on error
    }
}
void bs_tangent_t2n_seq_f(int start_index, unsigned int count, float* out) {
    if (!out || count == 0) return;
    try {
        boost::math::tangent_t2n<float>(start_index, count, out);
    } catch (...) {
        // Leave output as-is on error
    }
}
void bs_tangent_t2n_seq_l(int start_index, unsigned int count, long double* out) {
    if (!out || count == 0) return;
    try {
        boost::math::tangent_t2n<long double>(start_index, count, out);
    } catch (...) {
        // Leave output as-is on error
    }
}

// Fibonacci numbers
unsigned long long bs_fibonacci_ull(unsigned long long n) {
    return bs_wrap<unsigned long long>([&] { return boost::math::fibonacci<unsigned long long>(n); });
}

unsigned int bs_prime(unsigned int n) {
    return bs_wrap<unsigned int>([&] { return boost::math::prime(n); });
}

// Factorial
double bs_factorial(unsigned int i) {
    return bs_wrap<double>([&] { return boost::math::factorial<double>(i); });
}
float bs_factorial_f(unsigned int i) {
    return bs_wrap<float>([&] { return boost::math::factorial<float>(i); });
}
long double bs_factorial_l(unsigned int i) {
    return bs_wrap<long double>([&] { return boost::math::factorial<long double>(i); });
}

double bs_rising_factorial(double x, unsigned int i) {
    return bs_wrap<double>([&] { return boost::math::rising_factorial<double>(x, i); });
}
float bs_rising_factorial_f(float x, unsigned int i) {
    return bs_wrap<float>([&] { return boost::math::rising_factorial<float>(x, i); });
}
long double bs_rising_factorial_l(long double x, unsigned int i) {
    return bs_wrap<long double>([&] { return boost::math::rising_factorial<long double>(x, i); });
}

double bs_binomial_coefficient(unsigned int n, unsigned int k) {
    return bs_wrap<double>([&] { return boost::math::binomial_coefficient<double>(n, k); });
}
float bs_binomial_coefficient_f(unsigned int n, unsigned int k) {
    return bs_wrap<float>([&] { return boost::math::binomial_coefficient<float>(n, k); });
}
long double bs_binomial_coefficient_l(unsigned int n, unsigned int k) {
    return bs_wrap<long double>([&] { return boost::math::binomial_coefficient<long double>(n, k); });
}

double bs_double_factorial(unsigned int i) {
    return bs_wrap<double>([&] { return boost::math::double_factorial<double>(i); });
}
float bs_double_factorial_f(unsigned int i) {
    return bs_wrap<float>([&] { return boost::math::double_factorial<float>(i); });
}
long double bs_double_factorial_l(unsigned int i) {
    return bs_wrap<long double>([&] { return boost::math::double_factorial<long double>(i); });
}

} // extern "C"
