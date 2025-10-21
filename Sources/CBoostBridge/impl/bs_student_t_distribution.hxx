//
//  Created by Volker Thieme 2025.
//  Students t distribution bridge (Boost.Math) — implementation TU include
//
#include "../internal/bs_internal.hpp"
#include "../include/distributions/bs_student_t_distribution.h"

#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/complement.hpp>

using boost::math::students_t_distribution;
using boost::math::complement;


extern "C" {

// Static utility: required degrees of freedom
float bs_student_t_find_degrees_of_freedom_f(float difference_from_mean, float alpha, float beta, float sd, float hint) {
    return bs_wrap<float>([&]{ return students_t_distribution<float>::find_degrees_of_freedom(difference_from_mean, alpha, beta, sd, hint); });
}
double bs_student_t_find_degrees_of_freedom_d(double difference_from_mean, double alpha, double beta, double sd, double hint) {
    return bs_wrap<double>([&]{ return students_t_distribution<double>::find_degrees_of_freedom(difference_from_mean, alpha, beta, sd, hint); });
}
long double bs_student_t_find_degrees_of_freedom_l(long double difference_from_mean, long double alpha, long double beta, long double sd, long double hint) {
    return bs_wrap<long double>([&]{ return students_t_distribution<long double>::find_degrees_of_freedom(difference_from_mean, alpha, beta, sd, hint); });
}


} // extern "C"
