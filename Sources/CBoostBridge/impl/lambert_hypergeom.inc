// Lambert W and Hypergeometric functions
#include <boost/math/special_functions/lambert_w.hpp>
#include <boost/math/special_functions/hypergeometric_1F0.hpp>
#include <boost/math/special_functions/hypergeometric_0F1.hpp>
#include <boost/math/special_functions/hypergeometric_2F0.hpp>
#include <boost/math/special_functions/hypergeometric_1F1.hpp>
#include <boost/math/special_functions/hypergeometric_pFq.hpp>
#include <vector>
#include "../internal/bs_internal.hpp"

extern "C" {

// Lambert W (real branches)
double bs_lambert_w0(double x)        { return bs_wrap<double>([&] { return boost::math::lambert_w0(x); }); }
double bs_lambert_wm1(double x)       { return bs_wrap<double>([&] { return boost::math::lambert_wm1(x); }); }

float bs_lambert_w0_f(float x)        { return bs_wrap<float>([&] { return boost::math::lambert_w0(x); }); }
float bs_lambert_wm1_f(float x)       { return bs_wrap<float>([&] { return boost::math::lambert_wm1(x); }); }

long double bs_lambert_w0_l(long double x)  { return bs_wrap<long double>([&] { return boost::math::lambert_w0(x); }); }
long double bs_lambert_wm1_l(long double x) { return bs_wrap<long double>([&] { return boost::math::lambert_wm1(x); }); }

// Hypergeometric functions
double bs_hypergeometric_1F0(double a, double z) { return bs_wrap<double>([&] { return boost::math::hypergeometric_1F0(a, z); }); }
float  bs_hypergeometric_1F0_f(float a, float z) { return bs_wrap<float>([&] { return boost::math::hypergeometric_1F0(a, z); }); }
long double bs_hypergeometric_1F0_l(long double a, long double z) { return bs_wrap<long double>([&] { return boost::math::hypergeometric_1F0(a, z); }); }

double bs_hypergeometric_0F1(double b, double z) { return bs_wrap<double>([&] { return boost::math::hypergeometric_0F1(b, z); }); }
float  bs_hypergeometric_0F1_f(float b, float z) { return bs_wrap<float>([&] { return boost::math::hypergeometric_0F1(b, z); }); }
long double bs_hypergeometric_0F1_l(long double b, long double z) { return bs_wrap<long double>([&] { return boost::math::hypergeometric_0F1(b, z); }); }

double bs_hypergeometric_2F0(double a, double b, double z) { return bs_wrap<double>([&] { return boost::math::hypergeometric_2F0(a, b, z); }); }
float  bs_hypergeometric_2F0_f(float a, float b, float z)  { return bs_wrap<float>([&] { return boost::math::hypergeometric_2F0(a, b, z); }); }
long double bs_hypergeometric_2F0_l(long double a, long double b, long double z) { return bs_wrap<long double>([&] { return boost::math::hypergeometric_2F0(a, b, z); }); }

double bs_hypergeometric_1F1(double a, double b, double z) { return bs_wrap<double>([&] { return boost::math::hypergeometric_1F1(a, b, z); }); }
float  bs_hypergeometric_1F1_f(float a, float b, float z)  { return bs_wrap<float>([&] { return boost::math::hypergeometric_1F1(a, b, z); }); }
long double bs_hypergeometric_1F1_l(long double a, long double b, long double z) { return bs_wrap<long double>([&] { return boost::math::hypergeometric_1F1(a, b, z); }); }

double bs_hypergeometric_pFq(const double* a, size_t p, const double* b, size_t q, double z) {
    return bs_wrap<double>([&] {
        std::vector<double> va(a, a + p);
        std::vector<double> vb(b, b + q);
        return boost::math::hypergeometric_pFq(va, vb, z);
    });
}
float bs_hypergeometric_pFq_f(const float* a, size_t p, const float* b, size_t q, float z) {
    return bs_wrap<float>([&] {
        std::vector<float> va(a, a + p);
        std::vector<float> vb(b, b + q);
        return boost::math::hypergeometric_pFq(va, vb, z);
    });
}
long double bs_hypergeometric_pFq_l(const long double* a, size_t p, const long double* b, size_t q, long double z) {
    return bs_wrap<long double>([&] {
        std::vector<long double> va(a, a + p);
        std::vector<long double> vb(b, b + q);
        return boost::math::hypergeometric_pFq(va, vb, z);
    });
}

} // extern "C"

