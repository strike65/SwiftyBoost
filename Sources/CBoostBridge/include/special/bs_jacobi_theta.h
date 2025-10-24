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

double bs_jacobi_theta1_d(double k, double phi);
float bs_jacobi_theta1_f(float k, float phi);
long double bs_jacobi_theta1_l(long double k, long double phi);

double bs_jacobi_theta1tau_d(double k, double tau);
float bs_jacobi_theta1tau_f(float k, float tau);
long double bs_jacobi_theta1tau_l(long double k, long double tau);

double bs_jacobi_theta2_d(double k, double phi);
float bs_jacobi_theta2_f(float k, float phi);
long double bs_jacobi_theta2_l(long double k, long double phi);

double bs_jacobi_theta2tau_d(double k, double tau);
float bs_jacobi_theta2tau_f(float k, float tau);
long double bs_jacobi_theta2tau_l(long double k, long double tau);

double bs_jacobi_theta3_d(double k, double phi);
float bs_jacobi_theta3_f(float k, float phi);
long double bs_jacobi_theta3_l(long double k, long double phi);

double bs_jacobi_theta3tau_d(double k, double tau);
float bs_jacobi_theta3tau_f(float k, float tau);
long double bs_jacobi_theta3tau_l(long double k, long double tau);


double bs_jacobi_theta4_d(double k, double phi);
float bs_jacobi_theta4_f(float k, float phi);
long double bs_jacobi_theta4_l(long double k, long double phi);

double bs_jacobi_theta4tau_d(double k, double tau);
float bs_jacobi_theta4tau_f(float k, float tau);
long double bs_jacobi_theta4tau_l(long double k, long double tau);



#ifdef __cplusplus
}
#endif
