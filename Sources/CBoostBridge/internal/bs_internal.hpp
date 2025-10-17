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

