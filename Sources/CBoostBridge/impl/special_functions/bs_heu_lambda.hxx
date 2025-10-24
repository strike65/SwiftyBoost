// Heuman Lambda function wrappers
#include <boost/math/special_functions/heuman_lambda.hpp>
#include "../internal/bs_internal.hpp"

#ifdef __cplusplus
extern "C" {
#endif

double bs_heuman_lambda_d(double k, double phi) { return bs_wrap<double>([&] { return boost::math::heuman_lambda(k, phi); }); }
float bs_heuman_lambda_f(float k, float phi) { return bs_wrap<float>([&] { return boost::math::heuman_lambda(k, phi); }); }
long double bs_heuman_lambda_l(long double k, long double phi) { return bs_wrap<long double>([&] { return boost::math::heuman_lambda(k, phi); }); }

#ifdef __cplusplus
}
#endif
