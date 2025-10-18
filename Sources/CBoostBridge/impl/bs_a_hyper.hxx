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

#include <boost/math/special_functions/acosh.hpp>
#include <boost/math/special_functions/asinh.hpp>
#include <boost/math/special_functions/atanh.hpp>
#ifdef __cplusplus
extern "C"
{
#endif
    double bs_acosh(const double x) {
        return boost::math::acosh(x);
    }
    double bs_asinh(const double x) {
        return boost::math::asinh(x);
    }
    double bs_atanh(const double x) {
        return boost::math::atanh(x);
    }
    float bs_acosh_f(const float x) {
        return boost::math::atanh(x);
    }
    float bs_asinh_f(const float x) {
        return boost::math::atanh(x);
    }
    float bs_atanh_f(const float x) {
        return boost::math::atanh(x);
    }
    long double bs_acosh_l(const long double x) {
        return boost::math::atanh(x);
    }
    long double bs_asinh_l(const long double x) {
        return boost::math::atanh(x);
    }
    long double bs_atanh_l(const long double x) {
        return boost::math::atanh(x);
    }
#ifdef __cplusplus
}
#endif
