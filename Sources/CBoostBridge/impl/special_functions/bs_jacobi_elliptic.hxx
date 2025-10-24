//
//  Created by VT on 24.10.25.
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
#include <boost/math/special_functions/jacobi_elliptic.hpp>
#include "../internal/bs_internal.hpp"

extern "C"
{

double bs_jacobi_elliptic_d(double k, double theta, double* pcn, double* pdn) {
    double sn = bs_wrap<double>([&] { return boost::math::jacobi_elliptic(k, theta, pcn, pdn); });
    return sn;
}
float bs_jacobi_elliptic_f(float k, float theta, float* pcn, float* pdn) {
    float sn = bs_wrap<float>([&] { return boost::math::jacobi_elliptic(k, theta, pcn, pdn); });
    return sn;
}
long double bs_jacobi_elliptic_l(long double k, long double theta, long double* pcn, long double* pdn) {
    long double sn = bs_wrap<long double>([&] { return boost::math::jacobi_elliptic(k, theta, pcn, pdn); });
    return sn;
}

double bs_jacobi_elliptic_cd_d(double k, double theta) {
    return bs_wrap<double>([&] { return boost::math::jacobi_cd(k, theta); });
}
float bs_jacobi_elliptic_cd_f(float k, float theta) {
    return bs_wrap<float>([&] { return boost::math::jacobi_cd(k, theta); });
}
long double bs_jacobi_elliptic_cd_l(long double k, long double theta) {
    return bs_wrap<long double>([&] { return boost::math::jacobi_cd(k, theta); });
}

double bs_jacobi_elliptic_cn_d(double k, double theta) {
    return bs_wrap<double>([&] { return boost::math::jacobi_cn(k, theta); });
}
float bs_jacobi_elliptic_cn_f(float k, float theta) {
    return bs_wrap<float>([&] { return boost::math::jacobi_cn(k, theta); });
}
long double bs_jacobi_elliptic_cn_l(long double k, long double theta) {
    return bs_wrap<long double>([&] { return boost::math::jacobi_cn(k, theta); });
}

double bs_jacobi_elliptic_cs_d(double k, double theta) {
    return bs_wrap<double>([&] { return boost::math::jacobi_cs(k, theta); });
}
float bs_jacobi_elliptic_cs_f(float k, float theta) {
    return bs_wrap<float>([&] { return boost::math::jacobi_cs(k, theta); });
}
long double bs_jacobi_elliptic_cs_l(long double k, long double theta) {
    return bs_wrap<long double>([&] { return boost::math::jacobi_cs(k, theta); });
}

double bs_jacobi_elliptic_dc_d(double k, double theta) {
    return bs_wrap<double>([&] { return boost::math::jacobi_dc(k, theta); });
}
float bs_jacobi_elliptic_dc_f(float k, float theta) {
    return bs_wrap<float>([&] { return boost::math::jacobi_dc(k, theta); });
}
long double bs_jacobi_elliptic_dc_l(long double k, long double theta) {
    return bs_wrap<long double>([&] { return boost::math::jacobi_dc(k, theta); });
}

double bs_jacobi_elliptic_dn_d(double k, double theta) {
    return bs_wrap<double>([&] { return boost::math::jacobi_dn(k, theta); });
}
float bs_jacobi_elliptic_dn_f(float k, float theta) {
    return bs_wrap<float>([&] { return boost::math::jacobi_dn(k, theta); });
}
long double bs_jacobi_elliptic_dn_l(long double k, long double theta) {
    return bs_wrap<long double>([&] { return boost::math::jacobi_dn(k, theta); });
}

double bs_jacobi_elliptic_ds_d(double k, double theta) {
    return bs_wrap<double>([&] { return boost::math::jacobi_ds(k, theta); });
}
float bs_jacobi_elliptic_ds_f(float k, float theta) {
    return bs_wrap<float>([&] { return boost::math::jacobi_ds(k, theta); });
}
long double bs_jacobi_elliptic_ds_l(long double k, long double theta) {
    return bs_wrap<long double>([&] { return boost::math::jacobi_ds(k, theta); });
}

double bs_jacobi_elliptic_nc_d(double k, double theta) {
    return bs_wrap<double>([&] { return boost::math::jacobi_nc(k, theta); });
}
float bs_jacobi_elliptic_nc_f(float k, float theta) {
    return bs_wrap<float>([&] { return boost::math::jacobi_nc(k, theta); });
}
long double bs_jacobi_elliptic_nc_l(long double k, long double theta) {
    return bs_wrap<long double>([&] { return boost::math::jacobi_nc(k, theta); });
}

double bs_jacobi_elliptic_nd_d(double k, double theta) {
    return bs_wrap<double>([&] { return boost::math::jacobi_nd(k, theta); });
}
float bs_jacobi_elliptic_nd_f(float k, float theta) {
    return bs_wrap<float>([&] { return boost::math::jacobi_nd(k, theta); });
}
long double bs_jacobi_elliptic_nd_l(long double k, long double theta) {
    return bs_wrap<long double>([&] { return boost::math::jacobi_nd(k, theta); });
}

double bs_jacobi_elliptic_ns_d(double k, double theta) {
    return bs_wrap<double>([&] { return boost::math::jacobi_ns(k, theta); });
}
float bs_jacobi_elliptic_ns_f(float k, float theta) {
    return bs_wrap<float>([&] { return boost::math::jacobi_ns(k, theta); });
}
long double bs_jacobi_elliptic_ns_l(long double k, long double theta) {
    return bs_wrap<long double>([&] { return boost::math::jacobi_ns(k, theta); });
}

double bs_jacobi_elliptic_sc_d(double k, double theta) {
    return bs_wrap<double>([&] { return boost::math::jacobi_sc(k, theta); });
}
float bs_jacobi_elliptic_sc_f(float k, float theta) {
    return bs_wrap<float>([&] { return boost::math::jacobi_sc(k, theta); });
}
long double bs_jacobi_elliptic_sc_l(long double k, long double theta) {
    return bs_wrap<long double>([&] { return boost::math::jacobi_sc(k, theta); });
}

double bs_jacobi_elliptic_sd_d(double k, double theta) {
    return bs_wrap<double>([&] { return boost::math::jacobi_sd(k, theta); });
}
float bs_jacobi_elliptic_sd_f(float k, float theta) {
    return bs_wrap<float>([&] { return boost::math::jacobi_sd(k, theta); });
}
long double bs_jacobi_elliptic_sd_l(long double k, long double theta) {
    return bs_wrap<long double>([&] { return boost::math::jacobi_sd(k, theta); });
}

double bs_jacobi_elliptic_sn_d(double k, double theta) {
    return bs_wrap<double>([&] { return boost::math::jacobi_sn(k, theta); });
}
float bs_jacobi_elliptic_sn_f(float k, float theta) {
    return bs_wrap<float>([&] { return boost::math::jacobi_sn(k, theta); });
}
long double bs_jacobi_elliptic_sn_l(long double k, long double theta) {
    return bs_wrap<long double>([&] { return boost::math::jacobi_sn(k, theta); });
}

}
