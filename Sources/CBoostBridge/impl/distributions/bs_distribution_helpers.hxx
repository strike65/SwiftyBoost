//
//  Created by Volker Thieme 2025.
//  Students t distribution bridge (Boost.Math) — implementation TU include
//
#include "../internal/bs_internal.hpp"
#include "../include/distributions/bs_distribution_helpers.h"

#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/beta.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/complement.hpp>

using boost::math::students_t_distribution;
using boost::math::beta_distribution;
using boost::math::complement;
using boost::math::chi_squared_distribution;

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

// beta helpers
float bs_beta_find_alpha_f(float mean, float variance) {
    return bs_wrap<float>([&]{ return beta_distribution<float>::find_alpha(mean, variance); });
}
double bs_beta_find_alpha_d(double mean, double variance) {
    return bs_wrap<double>([&]{ return beta_distribution<double>::find_alpha(mean, variance); });
}
long double bs_beta_find_alpha_l(long double mean, long double variance) {
    return bs_wrap<long double>([&]{ return beta_distribution<long double>::find_alpha(mean, variance); });
}

float bs_beta_find_beta_f(float mean, float variance) {
    return bs_wrap<float>([&]{ return beta_distribution<float>::find_beta(mean, variance); });
}
double bs_beta_find_beta_d(double mean, double variance) {
    return bs_wrap<double>([&]{ return beta_distribution<double>::find_beta(mean, variance); });
}
long double bs_beta_find_beta_l(long double mean, long double variance) {
    return bs_wrap<long double>([&]{ return beta_distribution<long double>::find_beta(mean, variance); });
}
float bs_beta_find_alpha_from_beta_f(float beta, float x, float probability) {
    return bs_wrap<float>([&]{ return beta_distribution<float>::find_alpha(beta, x, probability); });
}
double bs_beta_find_alpha_from_beta_d(double beta, double x, double probability) {
    return bs_wrap<double>([&]{ return beta_distribution<double>::find_alpha(beta, x, probability); });
}
long double bs_beta_find_alpha_from_beta_l(long double beta, long double x, long double probability) {
    return bs_wrap<long double>([&]{ return beta_distribution<long double>::find_alpha(beta, x, probability); });
}

float bs_beta_find_beta_from_alpha_f(float alpha, float x, float probability) {
    return bs_wrap<float>([&]{ return beta_distribution<float>::find_alpha(alpha, x, probability); });
}
double bs_beta_find_beta_from_alpha_d(double alpha, double x, double probability) {
    return bs_wrap<double>([&]{ return beta_distribution<double>::find_alpha(alpha, x, probability); });
}
long double bs_beta_find_beta_from_alpha_l(long double alpha, long double x, long double probability) {
    return bs_wrap<long double>([&]{ return beta_distribution<long double>::find_alpha(alpha, x, probability); });
}
float bs_chisquare_find_degreesOfFreedom_f(float difference_from_variance, float alpha, float beta, float variance, float hint) {
    return bs_wrap<float>([&]{ return chi_squared_distribution<float>::find_degrees_of_freedom(difference_from_variance, alpha, beta, variance, hint); });
}

double bs_chisquare_find_degreesOfFreedom_d(double difference_from_variance, double alpha, double beta, double variance, double hint) {
    return bs_wrap<double>([&]{ return chi_squared_distribution<double>::find_degrees_of_freedom(difference_from_variance, alpha, beta, variance, hint); });
}
   
long double bs_chisquare_find_degreesOfFreedom_l(long double difference_from_variance, long double alpha, long double beta, long double variance, long double hint) {
    return bs_wrap<long double>([&]{ return chi_squared_distribution<long double>::find_degrees_of_freedom(difference_from_variance, alpha, beta, variance, hint); });
}



} // extern "C"
