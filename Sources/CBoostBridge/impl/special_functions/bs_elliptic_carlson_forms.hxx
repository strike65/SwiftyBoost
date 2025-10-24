// Elliptic integrals (Carlson symmetric forms)
#include <boost/math/special_functions/ellint_rc.hpp>
#include <boost/math/special_functions/ellint_rf.hpp>
#include <boost/math/special_functions/ellint_rd.hpp>
#include <boost/math/special_functions/ellint_rj.hpp>
#include <boost/math/special_functions/ellint_rg.hpp>
#include "../internal/bs_internal.hpp"

#ifdef __cplusplus
extern "C" {
#endif

double bs_ellint_rc_d(double x, double y)               { return bs_wrap<double>([&] { return boost::math::ellint_rc(x, y); }); }
double bs_ellint_rf_d(double x, double y, double z)     { return bs_wrap<double>([&] { return boost::math::ellint_rf(x, y, z); }); }
double bs_ellint_rd_d(double x, double y, double z)     { return bs_wrap<double>([&] { return boost::math::ellint_rd(x, y, z); }); }
double bs_ellint_rj_d(double x, double y, double z, double p) { return bs_wrap<double>([&] { return boost::math::ellint_rj(x, y, z, p); }); }
double bs_ellint_rg_d(double x, double y, double z)     { return bs_wrap<double>([&] { return boost::math::ellint_rg(x, y, z); }); }

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

#ifdef __cplusplus
}
#endif

