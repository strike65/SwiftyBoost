// Bessel (cylindrical and spherical) wrappers
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/bessel_prime.hpp>
#include "../internal/bs_internal.hpp"

#ifdef __cplusplus
extern "C" {
#endif

// Bessel (cylindrical, real)
double bs_cyl_bessel_j_d(double v, double x) { return bs_wrap<double>([&] { return boost::math::cyl_bessel_j(v, x); }); }
double bs_cyl_neumann_d(double v, double x)  { return bs_wrap<double>([&] { return boost::math::cyl_neumann(v, x); }); }
double bs_cyl_bessel_i_d(double v, double x) { return bs_wrap<double>([&] { return boost::math::cyl_bessel_i(v, x); }); }
double bs_cyl_bessel_k_d(double v, double x) { return bs_wrap<double>([&] { return boost::math::cyl_bessel_k(v, x); }); }
double bs_cyl_bessel_j_zero_d(double v, int m) { return bs_wrap<double>([&] { return boost::math::cyl_bessel_j_zero<double>(v, m); }); }
void bs_cyl_bessel_j_zeros_d(double v, int start_index, unsigned int number_of_zeros, double* out) {
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
double bs_sph_bessel_d(unsigned int n, double x) { return boost::math::sph_bessel(n, x); }
float bs_sph_bessel_f(unsigned int n, float x) { return boost::math::sph_bessel(n, x); }
long double bs_sph_bessel_l(unsigned int n, long double x) { return boost::math::sph_bessel(n, x); }
double bs_sph_neumann_d(unsigned int n, double x) { return boost::math::sph_neumann(n, x); }
float bs_sph_neumann_f(unsigned int n, float x) { return boost::math::sph_neumann(n, x); }
long double bs_sph_neumann_prime_l(unsigned int n, long double x); // fwd decl to keep order warnings away
long double bs_sph_neumann_l(unsigned int n, long double x) { return boost::math::sph_neumann(n, x); }

// Derivatives/primes
double bs_cyl_bessel_j_prime_d(double v, double x) { return bs_wrap<double>([&] { return boost::math::cyl_bessel_j_prime(v, x); }); }
double bs_cyl_bessel_i_prime_d(double v, double x) { return bs_wrap<double>([&] { return boost::math::cyl_bessel_i_prime(v, x); }); }
double bs_cyl_bessel_k_prime_d(double v, double x) { return bs_wrap<double>([&] { return boost::math::cyl_bessel_k_prime(v, x); }); }
double bs_sph_bessel_prime_d(unsigned int n, double x) { return bs_wrap<double>([&] { return boost::math::sph_bessel_prime(n, x); }); }
double bs_sph_neumann_prime_d(unsigned int n, double x) { return bs_wrap<double>([&] { return boost::math::sph_neumann_prime(n, x); }); }

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

#ifdef __cplusplus
}
#endif

