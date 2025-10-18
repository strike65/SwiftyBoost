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

// Prefer the generic header that provides real and complex sinc/sinhc
// Prefer the standard Boost.Math header sinc.hpp. For Boost versions that
// provide complex overloads in a separate header, include it conditionally.
#include <boost/math/special_functions/sinc.hpp>
#include <boost/math/special_functions/sinhc.hpp>
#include "../internal/bs_internal.hpp"
#include "../include/bs_complex.h"
#include <cmath>
#ifdef __cplusplus
extern "C" {
#endif
double bs_sinc_pi(double x) {
    return bs_wrap<double>([&] {
        if (x == 0.0) return 1.0;
        return boost::math::sin_pi(x) / (boost::math::constants::pi<double>() * x);
    });
}
bs_complex_d bs_sincc_pi(bs_complex_d x) {
    auto z = to_std(x);
    if (z == std::complex<double>(0.0, 0.0)) {
        return bs_complex_d{1.0, 0.0};
    }
    const double pi = boost::math::constants::pi<double>();
    auto denom = pi * z;
    auto y = std::sin(denom) / denom;
    return from_std(y);
}

float bs_sinc_pi_f(float x) {
    return bs_wrap<float>([&] {
        if (x == 0.0f) return 1.0f;
        const double xd = static_cast<double>(x);
        return static_cast<float>(boost::math::sin_pi(xd) / (boost::math::constants::pi<double>() * xd));
    });
}
bs_complex_f bs_sincc_pi_f(bs_complex_f x) {
    auto z = to_std(x);
    if (z == std::complex<float>(0.0f, 0.0f)) {
        return bs_complex_f{1.0f, 0.0f};
    }
    const float pi = static_cast<float>(boost::math::constants::pi<double>());
    auto denom = pi * z;
    auto y = std::sin(denom) / denom;
    return from_std(y);
}

long double bs_sinc_pi_l(long double x) {
    return bs_wrap<long double>([&] {
        if (x == static_cast<long double>(0.0)) return static_cast<long double>(1.0);
        return boost::math::sin_pi(x) / (boost::math::constants::pi<long double>() * x);
    });
}
bs_complex_l bs_sincc_pi_l(bs_complex_l x) {
    auto z = to_std(x);
    if (z == std::complex<long double>(0.0L, 0.0L)) {
        return bs_complex_l{static_cast<long double>(1.0), static_cast<long double>(0.0)};
    }
    const long double pi = boost::math::constants::pi<long double>();
    auto denom = pi * z;
    auto y = std::sin(denom) / denom;
    return from_std(y);
}
double bs_sinhc_pi(double x) {
    return bs_wrap<double>([&] {
        if (x == 0.0) return 1.0;
        const double p = boost::math::constants::pi<double>();
        return std::sinh(p * x) / (p * x);
    });
}
bs_complex_d bs_sinhcc_pi(bs_complex_d x) {
    auto y = from_std(boost::math::sinhc_pi(to_std(x)));
    bs_complex_d r { static_cast<double>(y.re), static_cast<double>(y.im) };
    return r;
}

float bs_sinhc_pi_f(float x) {
    return bs_wrap<float>([&] {
        if (x == 0.0f) return 1.0f;
        const double p = boost::math::constants::pi<double>();
        const double xd = static_cast<double>(x);
        return static_cast<float>(std::sinh(p * xd) / (p * xd));
    });
}
bs_complex_f bs_sinhcc_pi_f(bs_complex_f x) {
    auto y = from_std(boost::math::sinhc_pi(to_std(x)));
    bs_complex_f r { static_cast<float>(y.re), static_cast<float>(y.im) };
    return r;
}

long double bs_sinhc_pi_l(long double x) {
    return bs_wrap<long double>([&] {
        if (x == static_cast<long double>(0.0)) return static_cast<long double>(1.0);
        const long double p = boost::math::constants::pi<long double>();
        return std::sinh(p * x) / (p * x);
    });
}
bs_complex_l bs_sinhcc_pi_l(bs_complex_l x) {
    auto y = from_std(boost::math::sinhc_pi(to_std(x)));
    bs_complex_l r { static_cast<long double>(y.re), static_cast<long double>(y.im) };
    return r;
}

#ifdef __cplusplus
}
#endif
