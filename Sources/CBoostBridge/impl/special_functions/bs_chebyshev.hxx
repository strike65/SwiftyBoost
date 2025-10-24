// Chebyshev polynomials and series
#include <boost/math/special_functions/chebyshev.hpp>
#include "../internal/bs_internal.hpp"

extern "C" {

double bs_chebyshev_T_d(unsigned int n, double x) { return bs_wrap<double>([&] { return boost::math::chebyshev_t(n, x); }); }
double bs_chebyshev_U_d(unsigned int n, double x) { return bs_wrap<double>([&] { return boost::math::chebyshev_u(n, x); }); }
float bs_chebyshev_T_f(unsigned int n, float x) { return bs_wrap<float>([&] { return boost::math::chebyshev_t(n, x); }); }
float bs_chebyshev_U_f(unsigned int n, float x) { return bs_wrap<float>([&] { return boost::math::chebyshev_u(n, x); }); }
long double bs_chebyshev_T_l(unsigned int n, long double x) { return bs_wrap<long double>([&] { return boost::math::chebyshev_t(n, x); }); }
long double bs_chebyshev_U_l(unsigned int n, long double x) { return bs_wrap<long double>([&] { return boost::math::chebyshev_u(n, x); }); }

// Chebyshev series evaluation via Clenshaw (first kind)
double bs_chebyshev_clenshaw_d(const double* c, size_t count, double x) { return bs_wrap<double>([&] { return boost::math::chebyshev_clenshaw_recurrence(c, count, x);} ); }
float bs_chebyshev_clenshaw_f(const float* c, size_t count, float x) { return bs_wrap<float>([&] { return boost::math::chebyshev_clenshaw_recurrence(c, count, x);} ); }
long double bs_chebyshev_clenshaw_l(const long double* c, size_t count, long double x) { return bs_wrap<long double>([&] { return boost::math::chebyshev_clenshaw_recurrence(c, count, x);} ); }

} // extern "C"

