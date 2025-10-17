#pragma once

#ifdef __cplusplus
extern "C" {
#endif

// POD layout compatible with Swift's import.
typedef struct {
    float re, im;
} bs_complex_f;
typedef struct {
    double re, im;
} bs_complex_d;
typedef struct {
    long double re, im;
} bs_complex_l;

// Elementary complex arithmetic
bs_complex_d bs_cadd(bs_complex_d a, bs_complex_d b);
bs_complex_d bs_csub(bs_complex_d a, bs_complex_d b);
bs_complex_d bs_cmul(bs_complex_d a, bs_complex_d b);
bs_complex_d bs_cdiv(bs_complex_d a, bs_complex_d b);

bs_complex_f bs_cadd_f(bs_complex_f a, bs_complex_f b);
bs_complex_f bs_csub_f(bs_complex_f a, bs_complex_f b);
bs_complex_f bs_cmul_f(bs_complex_f a, bs_complex_f b);
bs_complex_f bs_cdiv_f(bs_complex_f a, bs_complex_f b);

bs_complex_l bs_cadd_l(bs_complex_l a, bs_complex_l b);
bs_complex_l bs_csub_l(bs_complex_l a, bs_complex_l b);
bs_complex_l bs_cmul_l(bs_complex_l a, bs_complex_l b);
bs_complex_l bs_cdiv_l(bs_complex_l a, bs_complex_l b);

// Elementary complex functions
bs_complex_d bs_cexp(bs_complex_d z);
bs_complex_d bs_clog(bs_complex_d z);
bs_complex_d bs_csqrt(bs_complex_d z);
bs_complex_d bs_csin(bs_complex_d z);
bs_complex_d bs_ccos(bs_complex_d z);
bs_complex_d bs_ctan(bs_complex_d z);
bs_complex_d bs_csinh(bs_complex_d z);
bs_complex_d bs_ccosh(bs_complex_d z);
bs_complex_d bs_ctanh(bs_complex_d z);
bs_complex_d bs_catan(bs_complex_d z);

bs_complex_f bs_cexp_f(bs_complex_f z);
bs_complex_f bs_clog_f(bs_complex_f z);
bs_complex_f bs_csqrt_f(bs_complex_f z);
bs_complex_f bs_csin_f(bs_complex_f z);
bs_complex_f bs_ccos_f(bs_complex_f z);
bs_complex_f bs_ctan_f(bs_complex_f z);
bs_complex_f bs_csinh_f(bs_complex_f z);
bs_complex_f bs_ccosh_f(bs_complex_f z);
bs_complex_f bs_ctanh_f(bs_complex_f z);
bs_complex_f bs_catan_f(bs_complex_f z);

bs_complex_l bs_cexp_l(bs_complex_l z);
bs_complex_l bs_clog_l(bs_complex_l z);
bs_complex_l bs_csqrt_l(bs_complex_l z);
bs_complex_l bs_csin_l(bs_complex_l z);
bs_complex_l bs_ccos_l(bs_complex_l z);
bs_complex_l bs_ctan_l(bs_complex_l z);
bs_complex_l bs_csinh_l(bs_complex_l z);
bs_complex_l bs_ccosh_l(bs_complex_l z);
bs_complex_l bs_ctanh_l(bs_complex_l z);
bs_complex_l bs_catan_l(bs_complex_l z);

#ifdef __cplusplus
}
#endif

