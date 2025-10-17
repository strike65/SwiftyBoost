// Bessel, Legendre polynomials, and Elliptic integrals
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/bessel_prime.hpp>
#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include <boost/math/special_functions/ellint_3.hpp>
#include <boost/math/special_functions/ellint_rc.hpp>
#include <boost/math/special_functions/ellint_rf.hpp>
#include <boost/math/special_functions/ellint_rd.hpp>
#include <boost/math/special_functions/ellint_rj.hpp>
#include <boost/math/special_functions/ellint_rg.hpp>
#include <algorithm>
#include <vector>
#include "../internal/bs_internal.hpp"

extern "C" {

// Bessel (cylindrical, real)
double bs_cyl_bessel_j(double v, double x) { return bs_wrap<double>([&] { return boost::math::cyl_bessel_j(v, x); }); }
double bs_cyl_neumann(double v, double x)  { return bs_wrap<double>([&] { return boost::math::cyl_neumann(v, x); }); }
double bs_cyl_bessel_i(double v, double x) { return bs_wrap<double>([&] { return boost::math::cyl_bessel_i(v, x); }); }
double bs_cyl_bessel_k(double v, double x) { return bs_wrap<double>([&] { return boost::math::cyl_bessel_k(v, x); }); }
double bs_cyl_bessel_j_zero(double v, int m) { return bs_wrap<double>([&] { return boost::math::cyl_bessel_j_zero<double>(v, m); }); }
void bs_cyl_bessel_j_zeros(double v, int start_index, unsigned int number_of_zeros, double* out) {
    if (!out || number_of_zeros == 0) return;
    boost::math::cyl_bessel_j_zero(v, start_index, number_of_zeros, out);
}

float bs_cyl_bessel_j_f(float v, float x) { return bs_wrap<float>([&] { return boost::math::cyl_bessel_j(v, x); }); }
float bs_cyl_neumann_f(float v, float x)  { return bs_wrap<float>([&] { return boost::math::cyl_neumann(v, x); }); }
float bs_cyl_bessel_i_f(float v, float x) { return bs_wrap<float>([&] { return boost::math::cyl_bessel_i(v, x); }); }
float bs_cyl_bessel_k_f(float v, float x) { return bs_wrap<float>([&] { return boost::math::cyl_bessel_k(v, x); }); }
float bs_cyl_bessel_j_zero_f(float v, int m) { return bs_wrap<float>([&] { return boost::math::cyl_bessel_j_zero<float>(v, m); }); }
void bs_cyl_bessel_j_zeros_f(float v, int start_index, unsigned int number_of_zeros, float* out) {
    if (!out || number_of_zeros == 0) return;
    boost::math::cyl_bessel_j_zero(v, start_index, number_of_zeros, out);
}

long double bs_cyl_bessel_j_l(long double v, long double x) { return bs_wrap<long double>([&] { return boost::math::cyl_bessel_j(v, x); }); }
long double bs_cyl_neumann_l(long double v, long double x)  { return bs_wrap<long double>([&] { return boost::math::cyl_neumann(v, x); }); }
long double bs_cyl_bessel_i_l(long double v, long double x) { return bs_wrap<long double>([&] { return boost::math::cyl_bessel_i(v, x); }); }
long double bs_cyl_bessel_k_l(long double v, long double x) { return bs_wrap<long double>([&] { return boost::math::cyl_bessel_k(v, x); }); }
long double bs_cyl_bessel_j_zero_l(long double v, int m) { return bs_wrap<long double>([&] { return boost::math::cyl_bessel_j_zero<long double>(v, m); }); }
void bs_cyl_bessel_j_zeros_l(long double v, int start_index, unsigned int number_of_zeros, long double* out) {
    if (!out || number_of_zeros == 0) return;
    boost::math::cyl_bessel_j_zero(v, start_index, number_of_zeros, out);
}

// Spherical Bessel/Neumann
double bs_sph_bessel(unsigned int n, double x) { return boost::math::sph_bessel(n, x); }
float bs_sph_bessel_f(unsigned int n, float x) { return boost::math::sph_bessel(n, x); }
long double bs_sph_bessel_l(unsigned int n, long double x) { return boost::math::sph_bessel(n, x); }

double bs_sph_neumann(unsigned int n, double x) { return boost::math::sph_neumann(n, x); }
float bs_sph_neumann_f(unsigned int n, float x) { return boost::math::sph_neumann(n, x); }
long double bs_sph_neumann_l(unsigned int n, long double x) { return boost::math::sph_neumann(n, x); }

double bs_cyl_bessel_j_prime(double v, double x) { return bs_wrap<double>([&] { return boost::math::cyl_bessel_j_prime(v, x); }); }
double bs_cyl_bessel_i_prime(double v, double x) { return bs_wrap<double>([&] { return boost::math::cyl_bessel_i_prime(v, x); }); }
double bs_cyl_bessel_k_prime(double v, double x) { return bs_wrap<double>([&] { return boost::math::cyl_bessel_k_prime(v, x); }); }
double bs_sph_bessel_prime(unsigned int n, double x) { return bs_wrap<double>([&] { return boost::math::sph_bessel_prime(n, x); }); }
double bs_sph_neumann_prime(unsigned int n, double x) { return bs_wrap<double>([&] { return boost::math::sph_neumann_prime(n, x); }); }

float bs_cyl_bessel_j_prime_f(float v, float x) { return bs_wrap<float>([&] { return boost::math::cyl_bessel_j_prime(v, x); }); }
float bs_cyl_bessel_i_prime_f(float v, float x) { return bs_wrap<float>([&] { return boost::math::cyl_bessel_i_prime(v, x); }); }
float bs_cyl_bessel_k_prime_f(float v, float x) { return bs_wrap<float>([&] { return boost::math::cyl_bessel_k_prime(v, x); }); }
float bs_sph_bessel_prime_f(unsigned int n, float x) { return bs_wrap<float>([&] { return boost::math::sph_bessel_prime(n, x); }); }
float bs_sph_neumann_prime_f(unsigned int n, float x) { return bs_wrap<float>([&] { return boost::math::sph_neumann_prime(n, x); }); }

long double bs_cyl_bessel_j_prime_l(long double v, long double x) { return bs_wrap<long double>([&] { return boost::math::cyl_bessel_j_prime(v, x); }); }
long double bs_cyl_bessel_i_prime_l(long double v, long double x) { return bs_wrap<long double>([&] { return boost::math::cyl_bessel_i_prime(v, x); }); }
long double bs_cyl_bessel_k_prime_l(long double v, long double x) { return bs_wrap<long double>([&] { return boost::math::cyl_bessel_k_prime(v, x); }); }
long double bs_sph_bessel_prime_l(unsigned int n, long double x) { return bs_wrap<long double>([&] { return boost::math::sph_bessel_prime(n, x); }); }
long double bs_sph_neumann_prime_l(unsigned int n, long double x) { return bs_wrap<long double>([&] { return boost::math::sph_neumann_prime(n, x); }); }

// Legendre
double bs_legendre_p(int n, double x)            { return bs_wrap<double>([&] { return boost::math::legendre_p(n, x); }); }
double bs_assoc_legendre_p(int n, int m, double x) { return bs_wrap<double>([&] { return boost::math::legendre_p(n, m, x); }); }
double bs_legendre_p_prime(int n, double x)      { return bs_wrap<double>([&] { return boost::math::legendre_p_prime(n, x); }); }
void bs_legendre_p_zeros(int l, double* out) {
    if (!out) return;
    try {
        auto v = boost::math::legendre_p_zeros<double>(l);
        const std::size_t n = std::min<std::size_t>(static_cast<std::size_t>(l), v.size());
        std::copy(v.begin(), v.begin() + n, out);
    } catch (...) {
        // Leave output as-is on error
    }
}

float bs_legendre_p_f(int n, float x)                 { return bs_wrap<float>([&] { return boost::math::legendre_p(n, x); }); }
float bs_assoc_legendre_p_f(int n, int m, float x)    { return bs_wrap<float>([&] { return boost::math::legendre_p(n, m, x); }); }
float bs_legendre_p_prime_f(int n, float x)           { return bs_wrap<float>([&] { return boost::math::legendre_p_prime(n, x); }); }
void bs_legendre_p_zeros_f(int l, float* out) {
    if (!out) return;
    try {
        auto v = boost::math::legendre_p_zeros<float>(l);
        const std::size_t n = std::min<std::size_t>(static_cast<std::size_t>(l), v.size());
        std::copy(v.begin(), v.begin() + n, out);
    } catch (...) {
        // Leave output as-is on error
    }
}

long double bs_legendre_p_l(int n, long double x)              { return bs_wrap<long double>([&] { return boost::math::legendre_p(n, x); }); }
long double bs_assoc_legendre_p_l(int n, int m, long double x) { return bs_wrap<long double>([&] { return boost::math::legendre_p(n, m, x); }); }
long double bs_legendre_p_prime_l(int n, long double x)        { return bs_wrap<long double>([&] { return boost::math::legendre_p_prime(n, x); }); }
void bs_legendre_p_zeros_l(int l, long double* out) {
    if (!out) return;
    try {
        auto v = boost::math::legendre_p_zeros<long double>(l);
        const std::size_t n = std::min<std::size_t>(static_cast<std::size_t>(l), v.size());
        std::copy(v.begin(), v.begin() + n, out);
    } catch (...) {
        // Leave output as-is on error
    }
}

// Elliptic integrals (Legendre forms)
double bs_ellint_1_complete(double k)                 { return bs_wrap<double>([&] { return boost::math::ellint_1(k); }); }
double bs_ellint_1(double k, double phi)              { return bs_wrap<double>([&] { return boost::math::ellint_1(k, phi); }); }
double bs_ellint_2_complete(double k)                 { return bs_wrap<double>([&] { return boost::math::ellint_2(k); }); }
double bs_ellint_2(double k, double phi)              { return bs_wrap<double>([&] { return boost::math::ellint_2(k, phi); }); }
double bs_ellint_3(double k, double nu, double phi)   { return bs_wrap<double>([&] { return boost::math::ellint_3(k, nu, phi); }); }
double bs_ellint_3_complete(double k, double nu)      { return bs_wrap<double>([&] { return boost::math::ellint_3(k, nu); }); }

float bs_ellint_1_complete_f(float k)                 { return bs_wrap<float>([&] { return boost::math::ellint_1(k); }); }
float bs_ellint_1_f(float k, float phi)               { return bs_wrap<float>([&] { return boost::math::ellint_1(k, phi); }); }
float bs_ellint_2_complete_f(float k)                 { return bs_wrap<float>([&] { return boost::math::ellint_2(k); }); }
float bs_ellint_2_f(float k, float phi)               { return bs_wrap<float>([&] { return boost::math::ellint_2(k, phi); }); }
float bs_ellint_3_f(float k, float nu, float phi)     { return bs_wrap<float>([&] { return boost::math::ellint_3(k, nu, phi); }); }
float bs_ellint_3_complete_f(float k, float nu)       { return bs_wrap<float>([&] { return boost::math::ellint_3(k, nu); }); }

long double bs_ellint_1_complete_l(long double k)                 { return bs_wrap<long double>([&] { return boost::math::ellint_1(k); }); }
long double bs_ellint_1_l(long double k, long double phi)         { return bs_wrap<long double>([&] { return boost::math::ellint_1(k, phi); }); }
long double bs_ellint_2_complete_l(long double k)                 { return bs_wrap<long double>([&] { return boost::math::ellint_2(k); }); }
long double bs_ellint_2_l(long double k, long double phi)         { return bs_wrap<long double>([&] { return boost::math::ellint_2(k, phi); }); }
long double bs_ellint_3_l(long double k, long double nu,long double phi) { return bs_wrap<long double>([&] { return boost::math::ellint_3(k, nu, phi); }); }
long double bs_ellint_3_complete_l(long double k, long double nu)      { return bs_wrap<long double>([&] { return boost::math::ellint_3(k, nu); }); }

// Elliptic integrals (Carlson symmetric forms)
double bs_ellint_rc(double x, double y)               { return bs_wrap<double>([&] { return boost::math::ellint_rc(x, y); }); }
double bs_ellint_rf(double x, double y, double z)     { return bs_wrap<double>([&] { return boost::math::ellint_rf(x, y, z); }); }
double bs_ellint_rd(double x, double y, double z)     { return bs_wrap<double>([&] { return boost::math::ellint_rd(x, y, z); }); }
double bs_ellint_rj(double x, double y, double z, double p) { return bs_wrap<double>([&] { return boost::math::ellint_rj(x, y, z, p); }); }
double bs_ellint_rg(double x, double y, double z)     { return bs_wrap<double>([&] { return boost::math::ellint_rg(x, y, z); }); }

float bs_ellint_rc_f(float x, float y)                { return bs_wrap<float>([&] { return boost::math::ellint_rc(x, y); }); }
float bs_ellint_rf_f(float x, float y, float z)       { return bs_wrap<float>([&] { return boost::math::ellint_rf(x, y, z); }); }
float bs_ellint_rd_f(float x, float y, float z)       { return bs_wrap<float>([&] { return boost::math::ellint_rd(x, y, z); }); }
float bs_ellint_rj_f(float x, float y, float z, float p) { return bs_wrap<float>([&] { return boost::math::ellint_rj(x, y, z, p); }); }
float bs_ellint_rg_f(float x, float y, float z)       { return bs_wrap<float>([&] { return boost::math::ellint_rg(x, y, z); }); }

long double bs_ellint_rc_l(long double x, long double y)                { return bs_wrap<long double>([&] { return boost::math::ellint_rc(x, y); }); }
long double bs_ellint_rf_l(long double x, long double y, long double z) { return bs_wrap<long double>([&] { return boost::math::ellint_rf(x, y, z); }); }
long double bs_ellint_rd_l(long double x, long double y, long double z) { return bs_wrap<long double>([&] { return boost::math::ellint_rd(x, y, z); }); }
long double bs_ellint_rj_l(long double x, long double y, long double z, long double p) { return bs_wrap<long double>([&] { return boost::math::ellint_rj(x, y, z, p); }); }
long double bs_ellint_rg_l(long double x, long double y, long double z) { return bs_wrap<long double>([&] { return boost::math::ellint_rg(x, y, z); }); }

} // extern "C"

