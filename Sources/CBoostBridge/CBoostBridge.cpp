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

// Boost constants
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

// Gamma / Error
double bs_tgamma(double x)              { return boost::math::tgamma(x); }
double bs_lgamma(double x)              { return boost::math::lgamma(x); }
double bs_erf(double x)                 { return boost::math::erf(x); }
double bs_erfc(double x)                { return boost::math::erfc(x); }

// Beta family
double bs_beta(double a, double b)                     { return boost::math::beta(a, b); }
double bs_fullBeta(double a, double b, double x)       { return boost::math::beta(a, b, x); }
double bs_ibeta(double a, double b, double x)          { return boost::math::ibeta(a, b, x); }
double bs_ibetac(double a, double b, double x)         { return boost::math::ibetac(a, b, x); }
double bs_ibeta_inv(double a, double b, double p)      { return boost::math::ibeta_inv(a, b, p); }
double bs_ibetac_inv(double a, double b, double p)     { return boost::math::ibetac_inv(a, b, p); }
double bs_ibeta_inva(double b, double x, double p)     { return boost::math::ibeta_inva(b, x, p); }
double bs_ibeta_invb(double a, double x, double p)     { return boost::math::ibeta_invb(a, x, p); }

// Digamma / Polygamma / Zeta
double bs_digamma(double x)            { return boost::math::digamma(x); }
double bs_trigamma(double x)           { return boost::math::trigamma(x); }
double bs_polygamma(int n, double x)   { return boost::math::polygamma(n, x); }
double bs_riemann_zeta(double x)       { return boost::math::zeta(x); }

// Owen's T
double bs_owens_t(double h, double a)  { return boost::math::owens_t(h, a); }

// Exponential integrals and related
double bs_expint_En(int n, double x)   { return boost::math::expint(n, x); }
double bs_expm1(double x)              { return boost::math::expm1(x); }
double bs_log1p(double x)              { return boost::math::log1p(x); }
double bs_log1pmx(double x)            { return boost::math::log1pmx(x); }
double bs_powm1(double x, double y)    { return boost::math::powm1(x, y); }
double bs_cbrt(double x)               { return boost::math::cbrt(x); }

// Trig helpers
double bs_sin_pi(double x)             { return boost::math::sin_pi(x); }
double bs_cos_pi(double x)             { return boost::math::cos_pi(x); }

// Airy
double bs_airy_ai(double x)            { return boost::math::airy_ai(x); }
double bs_airy_bi(double x)            { return boost::math::airy_bi(x); }
double bs_airy_ai_prime(double x)      { return boost::math::airy_ai_prime(x); }
double bs_airy_bi_prime(double x)      { return boost::math::airy_bi_prime(x); }

// Bessel (cylindrical, real)
double bs_cyl_bessel_j(double v, double x) { return boost::math::cyl_bessel_j(v, x); }
double bs_cyl_neumann(double v, double x)  { return boost::math::cyl_neumann(v, x); }
double bs_cyl_bessel_i(double v, double x) { return boost::math::cyl_bessel_i(v, x); }
double bs_cyl_bessel_k(double v, double x) { return boost::math::cyl_bessel_k(v, x); }

// Legendre
double bs_legendre_p(int n, double x)                 { return boost::math::legendre_p(n, x); }
double bs_assoc_legendre_p(int n, int m, double x)    { return boost::math::legendre_p(n, m, x); }

// Elliptic integrals (Legendre forms)
double bs_ellint_1_complete(double k)                 { return boost::math::ellint_1(k); }
double bs_ellint_1(double k, double phi)              { return boost::math::ellint_1(k, phi); }
double bs_ellint_2_complete(double k)                 { return boost::math::ellint_2(k); }
double bs_ellint_2(double k, double phi)              { return boost::math::ellint_2(k, phi); }
double bs_ellint_3_complete(double k, double nu)      { return boost::math::ellint_3(k, nu); }
double bs_ellint_3(double k, double nu, double phi)   { return boost::math::ellint_3(k, nu, phi); }

// Elliptic integrals (Carlson symmetric forms)
double bs_ellint_rc(double x, double y)               { return boost::math::ellint_rc(x, y); }
double bs_ellint_rf(double x, double y, double z)     { return boost::math::ellint_rf(x, y, z); }
double bs_ellint_rd(double x, double y, double z)     { return boost::math::ellint_rd(x, y, z); }
double bs_ellint_rj(double x, double y, double z, double p) { return boost::math::ellint_rj(x, y, z, p); }
double bs_ellint_rg(double x, double y, double z)     { return boost::math::ellint_rg(x, y, z); }

// Lambert W (real branches)
double bs_lambert_w0(double x)        { return boost::math::lambert_w0(x); }
double bs_lambert_wm1(double x)       { return boost::math::lambert_wm1(x); }

} // extern "C"
