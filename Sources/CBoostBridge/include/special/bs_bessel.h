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
#pragma once

#ifdef __cplusplus
extern "C" {
#endif

// Bessel (cylindrical, real)
double bs_cyl_bessel_j_d(double v, double x);
double bs_cyl_neumann_d(double v, double x);
double bs_cyl_bessel_i_d(double v, double x);
double bs_cyl_bessel_k_d(double v, double x);
double bs_cyl_bessel_j_zero_d(double v, int m);
void bs_cyl_bessel_j_zeros_d(double v, int start_index,
                           unsigned int number_of_zeros, double *out);

float bs_cyl_bessel_j_f(float v, float x);
float bs_cyl_neumann_f(float v, float x);
float bs_cyl_bessel_i_f(float v, float x);
float bs_cyl_bessel_k_f(float v, float x);
float bs_cyl_bessel_j_zero_f(float v, int m);
void bs_cyl_bessel_j_zeros_f(float v, int start_index,
                             unsigned int number_of_zeros, float *out);

long double bs_cyl_bessel_j_l(long double v, long double x);
long double bs_cyl_neumann_l(long double v, long double x);
long double bs_cyl_bessel_i_l(long double v, long double x);
long double bs_cyl_bessel_k_l(long double v, long double x);
long double bs_cyl_bessel_j_zero_l(long double v, int m);
void bs_cyl_bessel_j_zeros_l(long double v, int start_index,
                           unsigned int number_of_zeros, long double *out);

// Spherical Bessel/Neumann (real)
double bs_sph_bessel_d(unsigned int n, double x);
float bs_sph_bessel_f(unsigned int n, float x);
long double bs_sph_bessel_l(unsigned int n, long double x);
double bs_sph_neumann_d(unsigned int n, double x);
float bs_sph_neumann_f(unsigned int n, float x);
long double bs_sph_neumann_l(unsigned int n, long double x);
double bs_cyl_bessel_j_prime_d(double v, double x);
double bs_cyl_bessel_i_prime_d(double v, double x);
double bs_cyl_bessel_k_prime_d(double v, double x);
double bs_sph_bessel_prime_d(unsigned int n, double x);
double bs_sph_neumann_prime_d(unsigned int n, double x);

float bs_cyl_bessel_j_prime_f(float v, float x);
float bs_cyl_bessel_i_prime_f(float v, float x);
float bs_cyl_bessel_k_prime_f(float v, float x);
float bs_sph_bessel_prime_f(unsigned int n, float x);
float bs_sph_neumann_prime_f(unsigned int n, float x);

long double bs_cyl_bessel_j_prime_l(long double v, long double x);
long double bs_cyl_bessel_i_prime_l(long double v, long double x);
long double bs_cyl_bessel_k_prime_l(long double v, long double x);
long double bs_sph_bessel_prime_l(unsigned int n, long double x);
long double bs_sph_neumann_prime_l(unsigned int n, long double x);

#ifdef __cplusplus
}
#endif

