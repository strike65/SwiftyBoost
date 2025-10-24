// Numbers: Bernoulli, Tangent, Fibonacci, Prime
#include <boost/math/special_functions/bernoulli.hpp>
#include <boost/math/special_functions/fibonacci.hpp>
#include <boost/math/special_functions/prime.hpp>
#include <algorithm>
#include "../internal/bs_internal.hpp"

extern "C" {

// Bernoulli Numbers
double bs_bernoulli_b2n(const int n) { return bs_wrap<double>([&] { return boost::math::bernoulli_b2n<double>(n); }); }
float bs_bernoulli_b2n_f(const int n) { return bs_wrap<float>([&] { return boost::math::bernoulli_b2n<float>(n); }); }
long double bs_bernoulli_b2n_l(const int n) { return bs_wrap<long double>([&] { return boost::math::bernoulli_b2n<long double>(n); }); }

// Tangent Numbers (scalar)
double bs_tangent_t2n(const int n) { return bs_wrap<double>([&] { return boost::math::tangent_t2n<double>(n); }); }
float bs_tangent_t2n_f(const int n) { return bs_wrap<float>([&] { return boost::math::tangent_t2n<float>(n); }); }
long double bs_tangent_t2n_l(const int n) { return bs_wrap<long double>([&] { return boost::math::tangent_t2n<long double>(n); }); }

// Tangent Numbers (bulk sequence)
void bs_tangent_t2n_seq(int start_index, unsigned int count, double* out) {
    if (!out || count == 0) return;
    try {
        for (unsigned int i = 0; i < count; ++i) out[i] = boost::math::tangent_t2n<double>(start_index + static_cast<int>(i));
    } catch (...) {}
}
void bs_tangent_t2n_seq_f(int start_index, unsigned int count, float* out) {
    if (!out || count == 0) return;
    try {
        for (unsigned int i = 0; i < count; ++i) out[i] = boost::math::tangent_t2n<float>(start_index + static_cast<int>(i));
    } catch (...) {}
}
void bs_tangent_t2n_seq_l(int start_index, unsigned int count, long double* out) {
    if (!out || count == 0) return;
    try {
        for (unsigned int i = 0; i < count; ++i) out[i] = boost::math::tangent_t2n<long double>(start_index + static_cast<int>(i));
    } catch (...) {}
}

// Fibonacci numbers
unsigned long long bs_fibonacci_ull(unsigned long long n) { return bs_wrap<unsigned long long>([&] { return boost::math::fibonacci<unsigned long long>(n); }); }

// Prime numbers up to 10000th
unsigned int bs_prime(unsigned int n) { return bs_wrap<unsigned int>([&] { return boost::math::prime(n); }); }

} // extern "C"

