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

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

double bs_jacobi_elliptic_d(double k, double theta, double* pcn, double* pdn);
float bs_jacobi_elliptic_f(float k, float theta, float* pcn, float* pdn);
long double bs_jacobi_elliptic_l(long double k, long double theta, long double* pcn, long double* pdn);

double bs_jacobi_elliptic_cd_d(double k, double theta);
float bs_jacobi_elliptic_cd_f(float k, float theta);
long double bs_jacobi_elliptic_cd_l(long double k, long double theta);

double bs_jacobi_elliptic_cn_d(double k, double theta);
float bs_jacobi_elliptic_cn_f(float k, float theta);
long double bs_jacobi_elliptic_cn_l(long double k, long double theta);

double bs_jacobi_elliptic_cs_d(double k, double theta);
float bs_jacobi_elliptic_cs_f(float k, float theta);
long double bs_jacobi_elliptic_cs_l(long double k, long double theta);

double bs_jacobi_elliptic_dc_d(double k, double theta);
float bs_jacobi_elliptic_dc_f(float k, float theta);
long double bs_jacobi_elliptic_dc_l(long double k, long double theta);

double bs_jacobi_elliptic_dn_d(double k, double theta);
float bs_jacobi_elliptic_dn_f(float k, float theta);
long double bs_jacobi_elliptic_dn_l(long double k, long double theta);

double bs_jacobi_elliptic_ds_d(double k, double theta);
float bs_jacobi_elliptic_ds_f(float k, float theta);
long double bs_jacobi_elliptic_ds_l(long double k, long double theta);

double bs_jacobi_elliptic_nc_d(double k, double theta);
float bs_jacobi_elliptic_nc_f(float k, float theta);
long double bs_jacobi_elliptic_nc_l(long double k, long double theta);

double bs_jacobi_elliptic_nd_d(double k, double theta);
float bs_jacobi_elliptic_nd_f(float k, float theta);
long double bs_jacobi_elliptic_nd_l(long double k, long double theta);

double bs_jacobi_elliptic_ns_d(double k, double theta);
float bs_jacobi_elliptic_ns_f(float k, float theta);
long double bs_jacobi_elliptic_ns_l(long double k, long double theta);

double bs_jacobi_elliptic_sc_d(double k, double theta);
float bs_jacobi_elliptic_sc_f(float k, float theta);
long double bs_jacobi_elliptic_sc_l(long double k, long double theta);

double bs_jacobi_elliptic_sd_d(double k, double theta);
float bs_jacobi_elliptic_sd_f(float k, float theta);
long double bs_jacobi_elliptic_sd_l(long double k, long double theta);

double bs_jacobi_elliptic_sn_d(double k, double theta);
float bs_jacobi_elliptic_sn_f(float k, float theta);
long double bs_jacobi_elliptic_sn_l(long double k, long double theta);




#ifdef __cplusplus
}
#endif

