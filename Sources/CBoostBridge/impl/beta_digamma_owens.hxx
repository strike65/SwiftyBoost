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
// Beta family; Digamma/Polygamma/Zeta; Owen's T

#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>
#include <boost/math/special_functions/polygamma.hpp>
#include <boost/math/special_functions/zeta.hpp>
#include <boost/math/special_functions/owens_t.hpp>
#include "../internal/bs_internal.hpp"

#ifdef __cplusplus
extern "C" {
#endif
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
    
#ifdef __cplusplus
}
#endif

