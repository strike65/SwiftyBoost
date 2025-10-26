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
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/geometric.hpp>
#include <boost/math/distributions/negative_binomial.hpp>
#include <boost/math/distributions/non_central_chi_squared.hpp>
using boost::math::students_t_distribution;
using boost::math::beta_distribution;
using boost::math::complement;
using boost::math::chi_squared_distribution;
using boost::math::binomial_distribution;
using boost::math::geometric_distribution;
using boost::math::non_central_chi_squared_distribution;
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

float bs_binomial_find_lower_bound_on_p_f(float trials,
                                          float successes,
                                          float alpha,
                                          bool jeffreys) {
    if (jeffreys) {
        return bs_wrap<float>([&]{ return binomial_distribution<float>::find_lower_bound_on_p(trials, successes, alpha, binomial_distribution<float>::jeffreys_prior_interval);});
    }
    else {
        return bs_wrap<float>([&]{ return binomial_distribution<float>::find_lower_bound_on_p(trials, successes, alpha, binomial_distribution<float>::clopper_pearson_exact_interval); });
    }
}

double bs_binomial_find_lower_bound_on_p_d(double trials, double successes, double alpha, bool jeffreys) {
    if (jeffreys) {
        return bs_wrap<double>([&]{ return binomial_distribution<double>::find_lower_bound_on_p(trials, successes, alpha, binomial_distribution<double>::jeffreys_prior_interval); });
    }
    else {
        return bs_wrap<double>([&]{ return binomial_distribution<double>::find_lower_bound_on_p(trials, successes, alpha, binomial_distribution<double>::clopper_pearson_exact_interval); });
    }
}

long double bs_binomial_find_lower_bound_on_p_l(long double trials, long double successes, long double alpha, bool jeffreys) {
    if (jeffreys) {
        return bs_wrap<long double>([&]{ return binomial_distribution<long double>::find_lower_bound_on_p(trials, successes, alpha, binomial_distribution<long double>::jeffreys_prior_interval); });
    }
    else {
        return bs_wrap<long double>([&]{ return binomial_distribution<long double>::find_lower_bound_on_p(trials, successes, alpha, binomial_distribution<long double>::clopper_pearson_exact_interval); });
    }
}

float bs_binomial_find_upper_bound_on_p_f(float trials,
                                          float successes,
                                          float alpha,
                                          bool jeffreys) {
    if (jeffreys) {
        return bs_wrap<float>([&]{ return binomial_distribution<float>::find_upper_bound_on_p(trials, successes, alpha, binomial_distribution<float>::jeffreys_prior_interval); });
    }
    else {
        return bs_wrap<float>([&]{ return binomial_distribution<float>::find_upper_bound_on_p(trials, successes, alpha, binomial_distribution<float>::clopper_pearson_exact_interval); });
    }
}

double bs_binomial_find_upper_bound_on_p_d(double trials, double successes, double alpha, bool jeffreys) {
    if (jeffreys) {
        return bs_wrap<double>([&]{ return binomial_distribution<double>::find_upper_bound_on_p(trials, successes, alpha, binomial_distribution<double>::jeffreys_prior_interval); });
    }
    else {
        return bs_wrap<double>([&]{ return binomial_distribution<double>::find_upper_bound_on_p(trials, successes, alpha, binomial_distribution<double>::clopper_pearson_exact_interval); });
    }
}

long double bs_binomial_find_upper_bound_on_p_l(long double trials, long double successes, long double alpha, bool jeffreys) {
    if (jeffreys) {
        return bs_wrap<long double>([&]{ return binomial_distribution<long double>::find_upper_bound_on_p(trials, successes, alpha, binomial_distribution<long double>::jeffreys_prior_interval); });
    }
    else {
        return bs_wrap<long double>([&]{ return binomial_distribution<long double>::find_upper_bound_on_p(trials, successes, alpha, binomial_distribution<long double>::clopper_pearson_exact_interval); });
    }
}

float bs_binomial_find_minimum_number_of_trials_f(float events, float success_fraction, float alpha) {
    return bs_wrap<float>([&]{ return binomial_distribution<float>::find_minimum_number_of_trials(events, success_fraction, alpha); });
}

double bs_binomial_find_minimum_number_of_trials_d(double events, double success_fraction, double alpha){
    return bs_wrap<double>([&]{ return binomial_distribution<double>::find_minimum_number_of_trials(events, success_fraction, alpha); });
}
long double bs_binomial_find_minimum_number_of_trials_l(long double events, long double success_fraction, long double alpha){
    return bs_wrap<long double>([&]{ return binomial_distribution<long double>::find_minimum_number_of_trials(events, success_fraction, alpha); });
}
float bs_binomial_find_maximum_number_of_trials_f(float events, float success_fraction, float alpha) {
    return bs_wrap<float>([&]{ return binomial_distribution<float>::find_maximum_number_of_trials(events, success_fraction, alpha); });
}
double bs_binomial_find_maximum_number_of_trials_d(double events, double success_fraction, double alpha) {
    return bs_wrap<double>([&]{ return binomial_distribution<double>::find_maximum_number_of_trials(events, success_fraction, alpha); });
}

long double bs_binomial_find_maximum_number_of_trials_l(long double events, long double success_fraction, long double alpha) {
    return bs_wrap<long double>([&]{ return binomial_distribution<long double>::find_maximum_number_of_trials(events, success_fraction, alpha); });
}

// geometric

float bs_geometric_find_lower_bound_on_p_f(float trials,
                                          float alpha) {
    return bs_wrap<float>([&]{ return geometric_distribution<float>::find_lower_bound_on_p(trials, alpha);});
}

double bs_geometric_find_lower_bound_on_p_d(double trials, double alpha) {
    return bs_wrap<double>([&]{ return geometric_distribution<double>::find_lower_bound_on_p(trials, alpha);});
}

long double bs_geometric_find_lower_bound_on_p_l(long double trials, long double alpha) {
    return bs_wrap<long double>([&]{ return geometric_distribution<long double>::find_lower_bound_on_p(trials, alpha);});
}

float bs_geometric_find_upper_bound_on_p_f(float trials,
                                          float alpha) {
    return bs_wrap<float>([&]{ return geometric_distribution<float>::find_upper_bound_on_p(trials, alpha);});
}

double bs_geometric_find_upper_bound_on_p_d(double trials, double alpha) {
    return bs_wrap<float>([&]{ return geometric_distribution<float>::find_upper_bound_on_p(trials, alpha);});
}

long double bs_geometric_find_upper_bound_on_p_l(long double trials, long double alpha) {
    return bs_wrap<long double>([&]{ return geometric_distribution<long double>::find_upper_bound_on_p(trials, alpha);});
}


float bs_geometric_find_minimum_number_of_trials_f(float failures, float success_fraction, float alpha) {
    return bs_wrap<float>([&]{ return geometric_distribution<float>::find_minimum_number_of_trials(failures, success_fraction, alpha); });
}

double bs_geometric_find_minimum_number_of_trials_d(double failures, double success_fraction, double alpha){
    return bs_wrap<double>([&]{ return geometric_distribution<double>::find_minimum_number_of_trials(failures, success_fraction, alpha); });
}
long double bs_geometric_find_minimum_number_of_trials_l(long double failures, long double success_fraction, long double alpha){
    return bs_wrap<long double>([&]{ return geometric_distribution<long double>::find_minimum_number_of_trials(failures, success_fraction, alpha); });
}
float bs_geometric_find_maximum_number_of_trials_f(float failures, float success_fraction, float alpha) {
    return bs_wrap<float>([&]{ return geometric_distribution<float>::find_maximum_number_of_trials(failures, success_fraction, alpha); });
}

double bs_geometric_find_maximum_number_of_trials_d(double failures, double success_fraction, double alpha){
    return bs_wrap<double>([&]{ return geometric_distribution<double>::find_maximum_number_of_trials(failures, success_fraction, alpha); });
}
long double bs_geometric_find_maximum_number_of_trials_l(long double failures, long double success_fraction, long double alpha){
    return bs_wrap<long double>([&]{ return geometric_distribution<long double>::find_maximum_number_of_trials(failures, success_fraction, alpha); });
}

float bs_negative_binomial_find_lower_bound_on_p_f(float trials,
                                          float successes,
                                          float alpha) {
        return bs_wrap<float>([&]{ return negative_binomial_distribution<float>::find_lower_bound_on_p(trials, successes, alpha);});
}

double bs_negative_binomial_find_lower_bound_on_p_d(double trials, double successes, double alpha) {
        return bs_wrap<double>([&]{ return negative_binomial_distribution<double>::find_lower_bound_on_p(trials, successes, alpha); });
}

long double bs_negative_binomial_find_lower_bound_on_p_l(long double trials, long double successes, long double alpha) {
        return bs_wrap<long double>([&]{ return negative_binomial_distribution<long double>::find_lower_bound_on_p(trials, successes, alpha); });
}

float bs_negative_binomial_find_upper_bound_on_p_f(float trials,
                                          float successes,
                                          float alpha) {
        return bs_wrap<float>([&]{ return negative_binomial_distribution<float>::find_upper_bound_on_p(trials, successes, alpha);});
}

double bs_negative_binomial_find_upper_bound_on_p_d(double trials, double successes, double alpha) {
        return bs_wrap<double>([&]{ return negative_binomial_distribution<double>::find_upper_bound_on_p(trials, successes, alpha); });
}

long double bs_negative_binomial_find_upper_bound_on_p_l(long double trials, long double successes, long double alpha) {
        return bs_wrap<long double>([&]{ return negative_binomial_distribution<long double>::find_upper_bound_on_p(trials, successes, alpha); });
}

float bs_negative_binomial_find_minimum_number_of_trials_f(float failures, float success_fraction, float alpha) {
    return bs_wrap<float>([&]{ return negative_binomial_distribution<float>::find_minimum_number_of_trials(failures, success_fraction, alpha); });
}

double bs_negative_binomial_find_minimum_number_of_trials_d(double failures, double success_fraction, double alpha){
    return bs_wrap<double>([&]{ return negative_binomial_distribution<double>::find_minimum_number_of_trials(failures, success_fraction, alpha); });
}
long double bs_negative_binomial_find_minimum_number_of_trials_l(long double failures, long double success_fraction, long double alpha){
    return bs_wrap<long double>([&]{ return negative_binomial_distribution<long double>::find_minimum_number_of_trials(failures, success_fraction, alpha); });
}

float bs_negative_binomial_find_maximum_number_of_trials_f(float failures, float success_fraction, float alpha) {
    return bs_wrap<float>([&]{ return negative_binomial_distribution<float>::find_maximum_number_of_trials(failures, success_fraction, alpha); });
}

double bs_negative_binomial_find_maximum_number_of_trials_d(double failures, double success_fraction, double alpha){
    return bs_wrap<double>([&]{ return negative_binomial_distribution<double>::find_maximum_number_of_trials(failures, success_fraction, alpha); });
}
long double bs_negative_binomial_find_maximum_number_of_trials_l(long double failures, long double success_fraction, long double alpha){
    return bs_wrap<long double>([&]{ return negative_binomial_distribution<long double>::find_maximum_number_of_trials(failures, success_fraction, alpha); });
}

float bs_non_central_chisquare_find_degreesOfFreedom_f(float lambda, float x, float p) {
    return bs_wrap<float>([&]{ return non_central_chi_squared_distribution<float>::find_degrees_of_freedom(lambda, x, p); });
}
double bs_non_central_chisquare_find_degreesOfFreedom_d(double lambda, double x, double p) {
    return bs_wrap<double>([&]{ return non_central_chi_squared_distribution<double>::find_degrees_of_freedom(lambda, x, p); });
}
float bs_non_central_chisquare_find_degreesOfFreedom_l(long double lambda, long double x, long double p) {
    return bs_wrap<long double>([&]{ return non_central_chi_squared_distribution<long double>::find_degrees_of_freedom(lambda, x, p); });
}

float bs_non_central_chisquare_find_non_centrality_f(float v, float x, float p) {
    return bs_wrap<float>([&]{ return non_central_chi_squared_distribution<float>::find_non_centrality(v, x, p); });
}
double bs_non_central_chisquare_find_non_centrality_d(double v, double x, double p) {
    return bs_wrap<double>([&]{ return non_central_chi_squared_distribution<double>::find_non_centrality(v, x, p); });
}
long double bs_non_central_chisquare_find_non_centrality_l(long double v, long double x, long double p) {
    return bs_wrap<long double>([&]{ return non_central_chi_squared_distribution<long double>::find_non_centrality(v, x, p); });
}


} // extern "C"
