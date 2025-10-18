//
//  Created by Volker Thieme 2025.
//  Copyright © 2025 Volker Thieme.
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
// Constants (Boost.Math constants)
#include <boost/math/constants/constants.hpp>
#include "../internal/bs_internal.hpp"

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

} // extern "C"

