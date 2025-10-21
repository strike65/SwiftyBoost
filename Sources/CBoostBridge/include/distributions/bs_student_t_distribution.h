//
//  Created by Volker Thieme 2025.
//  Students t distribution bridge (Boost.Math) — C ABI
//
#pragma once

#ifdef __cplusplus
extern "C" {
#endif

// Planning utility: find required degrees of freedom.
// Mirrors boost::math::students_t_distribution<>::find_degrees_of_freedom.
// Parameters:
//  - difference_from_mean: effect size (difference in means)
//  - alpha: significance level (Type I error)
//  - beta: Type II error (1 - power)
//  - sd: standard deviation
//  - hint: initial guess for df (<= 0 -> default 1)
float       bs_student_t_find_degrees_of_freedom_f(float difference_from_mean, float alpha, float beta, float sd, float hint);
double      bs_student_t_find_degrees_of_freedom_d(double difference_from_mean, double alpha, double beta, double sd, double hint);
long double bs_student_t_find_degrees_of_freedom_l(long double difference_from_mean, long double alpha, long double beta, long double sd, long double hint);

#ifdef __cplusplus
}
#endif
