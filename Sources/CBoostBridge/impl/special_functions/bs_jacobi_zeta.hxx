// Jacobi Zeta function wrapper
#include <boost/math/special_functions/jacobi_zeta.hpp>
#include "../internal/bs_internal.hpp"

extern "C" {

double bs_jacobi_zeta_d(double k, double phi) {
    return bs_wrap<double>([&] { return boost::math::jacobi_zeta(k, phi); });
}

float bs_jacobi_zeta_f(float k, float phi) {
    return bs_wrap<float>([&] { return boost::math::jacobi_zeta(k, phi); });
}

long double bs_jacobi_zeta_l(long double k, long double phi) {
    return bs_wrap<long double>([&] { return boost::math::jacobi_zeta(k, phi); });
}

} // extern "C"
