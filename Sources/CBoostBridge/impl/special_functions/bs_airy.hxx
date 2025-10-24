// Airy functions
#include <boost/math/special_functions/airy.hpp>
#include "../internal/bs_internal.hpp"

#ifdef __cplusplus
extern "C" {
#endif

double bs_airy_ai_d(double x)            { return bs_wrap<double>([&] { return boost::math::airy_ai(x); }); }
double bs_airy_bi_d(double x)            { return bs_wrap<double>([&] { return boost::math::airy_bi(x); }); }
double bs_airy_ai_prime_d(double x)      { return bs_wrap<double>([&] { return boost::math::airy_ai_prime(x); }); }
double bs_airy_bi_prime_d(double x)      { return bs_wrap<double>([&] { return boost::math::airy_bi_prime(x); }); }
double bs_airy_ai_zero_d(int m)          { return bs_wrap<double>([&] { return boost::math::airy_ai_zero<double>(m + 1); }); }
double bs_airy_bi_zero_d(int m)          { return bs_wrap<double>([&] { return boost::math::airy_bi_zero<double>(m + 1); }); }
void bs_airy_ai_zeros_d(int start_index, unsigned int number_of_zeros, double* out) {
    if (!out || number_of_zeros == 0) return;
    for (unsigned int i = 0; i < number_of_zeros; ++i) {
        const int m = start_index + static_cast<int>(i) + 1;
        out[i] = bs_wrap<double>([&] { return boost::math::airy_ai_zero<double>(m); });
    }
}
void bs_airy_bi_zeros_d(int start_index, unsigned int number_of_zeros, double* out) {
    if (!out || number_of_zeros == 0) return;
    for (unsigned int i = 0; i < number_of_zeros; ++i) {
        const int m = start_index + static_cast<int>(i) + 1;
        out[i] = bs_wrap<double>([&] { return boost::math::airy_bi_zero<double>(m); });
    }
}

float bs_airy_ai_f(float x)            { return bs_wrap<float>([&] { return boost::math::airy_ai(x); }); }
float bs_airy_bi_f(float x)            { return bs_wrap<float>([&] { return boost::math::airy_bi(x); }); }
float bs_airy_ai_prime_f(float x)      { return bs_wrap<float>([&] { return boost::math::airy_ai_prime(x); }); }
float bs_airy_bi_prime_f(float x)      { return bs_wrap<float>([&] { return boost::math::airy_bi_prime(x); }); }
float bs_airy_ai_zero_f(int m)          { return bs_wrap<float>([&] { return boost::math::airy_ai_zero<float>(m + 1); }); }
float bs_airy_bi_zero_f(int m)          { return bs_wrap<float>([&] { return boost::math::airy_bi_zero<float>(m + 1); }); }
void bs_airy_ai_zeros_f(int start_index, unsigned int number_of_zeros, float* out) {
    if (!out || number_of_zeros == 0) return;
    for (unsigned int i = 0; i < number_of_zeros; ++i) {
        const int m = start_index + static_cast<int>(i) + 1;
        out[i] = bs_wrap<float>([&] { return boost::math::airy_ai_zero<float>(m); });
    }
}
void bs_airy_bi_zeros_f(int start_index, unsigned int number_of_zeros, float* out) {
    if (!out || number_of_zeros == 0) return;
    for (unsigned int i = 0; i < number_of_zeros; ++i) {
        const int m = start_index + static_cast<int>(i) + 1;
        out[i] = bs_wrap<float>([&] { return boost::math::airy_bi_zero<float>(m); });
    }
}

long double bs_airy_ai_l(long double x)       { return bs_wrap<long double>([&] { return boost::math::airy_ai(x); }); }
long double bs_airy_bi_l(long double x)       { return bs_wrap<long double>([&] { return boost::math::airy_bi(x); }); }
long double bs_airy_ai_prime_l(long double x) { return bs_wrap<long double>([&] { return boost::math::airy_ai_prime(x); }); }
long double bs_airy_bi_prime_l(long double x) { return bs_wrap<long double>([&] { return boost::math::airy_bi_prime(x); }); }
long double bs_airy_ai_zero_l(int m)          { return bs_wrap<long double>([&] { return boost::math::airy_ai_zero<long double>(m + 1); }); }
long double bs_airy_bi_zero_l(int m)          { return bs_wrap<long double>([&] { return boost::math::airy_bi_zero<long double>(m + 1); }); }
void bs_airy_ai_zeros_l(int start_index, unsigned int number_of_zeros, long double* out) {
    if (!out || number_of_zeros == 0) return;
    for (unsigned int i = 0; i < number_of_zeros; ++i) {
        const int m = start_index + static_cast<int>(i) + 1;
        out[i] = bs_wrap<long double>([&] { return boost::math::airy_ai_zero<long double>(m); });
    }
}
void bs_airy_bi_zeros_l(int start_index, unsigned int number_of_zeros, long double* out) {
    if (!out || number_of_zeros == 0) return;
    for (unsigned int i = 0; i < number_of_zeros; ++i) {
        const int m = start_index + static_cast<int>(i) + 1;
        out[i] = bs_wrap<long double>([&] { return boost::math::airy_bi_zero<long double>(m); });
    }
}

#ifdef __cplusplus
}
#endif
