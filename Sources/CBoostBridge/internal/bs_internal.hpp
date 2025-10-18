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
// Internal helpers shared across CBoostBridge implementation files.
#pragma once

#include <limits>
#include <stdexcept>
#include <complex>

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

// Complex wrapper with same policy: on overflow => {+Inf,+Inf}, on other => {NaN,NaN}.
template <class T, class F>
static inline std::complex<T> bs_wrap_complex(F&& f) noexcept {
    try {
        return f();
    } catch (const std::overflow_error&) {
        const T inf = std::numeric_limits<T>::infinity();
        return std::complex<T>(inf, inf);
    } catch (const std::domain_error&) {
        const T nan = std::numeric_limits<T>::quiet_NaN();
        return std::complex<T>(nan, nan);
    } catch (...) {
        const T nan = std::numeric_limits<T>::quiet_NaN();
        return std::complex<T>(nan, nan);
    }
}

