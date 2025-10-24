// Owen's T wrappers
#include <boost/math/special_functions/owens_t.hpp>
#include "../internal/bs_internal.hpp"

#ifdef __cplusplus
extern "C" {
#endif

double bs_owens_t_d(double h, double a)           { return bs_wrap<double>([&] { return boost::math::owens_t(h, a); }); }
float  bs_owens_t_f(float h, float a)             { return bs_wrap<float>([&] { return boost::math::owens_t(h, a); }); }
long double bs_owens_t_l(long double h, long double a) { return bs_wrap<long double>([&] { return boost::math::owens_t(h, a); }); }

#ifdef __cplusplus
}
#endif

