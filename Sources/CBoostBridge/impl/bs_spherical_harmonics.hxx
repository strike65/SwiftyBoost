// Spherical harmonics
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include "../internal/bs_internal.hpp"
#include "complex.h"

extern "C" {

complex_d bs_spherical_harmonic_d(unsigned int n, int m, double theta, double phi) {
    auto z = boost::math::spherical_harmonic(n, m, theta, phi);
    complex_d r{ z.real(), z.imag() };
    return r;
}
complex_f bs_spherical_harmonic_f(unsigned int n, int m, float theta, float phi) {
    auto z = boost::math::spherical_harmonic(n, m, theta, phi);
    complex_f r{ static_cast<float>(z.real()), static_cast<float>(z.imag()) };
    return r;
}
complex_l bs_spherical_harmonic_l(unsigned int n, int m, long double theta, long double phi) {
    auto z = boost::math::spherical_harmonic(n, m, theta, phi);
    complex_l r{ static_cast<long double>(z.real()), static_cast<long double>(z.imag()) };
    return r;
}

} // extern "C"

