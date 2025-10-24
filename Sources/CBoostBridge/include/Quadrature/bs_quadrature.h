#ifndef BOOST_QUADRATURE_WRAPPER_HPP
#define BOOST_QUADRATURE_WRAPPER_HPP

#ifdef __cplusplus
extern "C" {
#endif

// Opaque Handles
typedef void* QuadratureHandle;

// Function pointers for integrands per precision
typedef float (*IntegrandFunctionF)(float x, void* context);
typedef double (*IntegrandFunctionD)(double x, void* context);
typedef long double (*IntegrandFunctionL)(long double x, void* context);

// Quadrature types
typedef enum {
    QUAD_GAUSS_LEGENDRE = 0,
    QUAD_GAUSS_HERMITE = 1,
    QUAD_GAUSS_LAGUERRE = 2,
    QUAD_GAUSS_JACOBI = 3,
    QUAD_GAUSS_KRONROD = 4,
    QUAD_TANH_SINH = 5,
    QUAD_SINH_SINH = 6,
    QUAD_EXP_SINH = 7
} QuadratureType;

typedef enum {
    QUAD_PRECISION_FLOAT = 0,
    QUAD_PRECISION_DOUBLE = 1,
    QUAD_PRECISION_LONG_DOUBLE = 2
} QuadraturePrecision;

// Result structures per precision
typedef struct {
    float result;
    float error;
    float l1_norm;
    int iterations;
    int function_calls;
    int converged; // 1 = true, 0 = false
} QuadratureResultF;

typedef struct {
    double result;
    double error;
    double l1_norm;
    int iterations;
    int function_calls;
    int converged; // 1 = true, 0 = false
} QuadratureResultD;

typedef struct {
    long double result;
    long double error;
    long double l1_norm;
    int iterations;
    int function_calls;
    int converged; // 1 = true, 0 = false
} QuadratureResultL;

// === Gauss-Legendre (Standard) ===
QuadratureHandle quad_gauss_create_f(int points);
QuadratureHandle quad_gauss_create_d(int points);
QuadratureHandle quad_gauss_create_l(int points);

// === Gauss-Kronrod ===
QuadratureHandle quad_gauss_kronrod_create_f(int points);
QuadratureHandle quad_gauss_kronrod_create_d(int points);
QuadratureHandle quad_gauss_kronrod_create_l(int points);

// === Tanh-Sinh (Double Exponential) ===
QuadratureHandle quad_tanh_sinh_create_f(void);
QuadratureHandle quad_tanh_sinh_create_d(void);
QuadratureHandle quad_tanh_sinh_create_l(void);

QuadratureHandle quad_tanh_sinh_create_with_params_f(int max_refinements, float tolerance);
QuadratureHandle quad_tanh_sinh_create_with_params_d(int max_refinements, double tolerance);
QuadratureHandle quad_tanh_sinh_create_with_params_l(int max_refinements, long double tolerance);

// === Sinh-Sinh ===
QuadratureHandle quad_sinh_sinh_create_f(void);
QuadratureHandle quad_sinh_sinh_create_d(void);
QuadratureHandle quad_sinh_sinh_create_l(void);

QuadratureHandle quad_sinh_sinh_create_with_params_f(int max_refinements, float tolerance);
QuadratureHandle quad_sinh_sinh_create_with_params_d(int max_refinements, double tolerance);
QuadratureHandle quad_sinh_sinh_create_with_params_l(int max_refinements, long double tolerance);

// === Exp-Sinh ===
QuadratureHandle quad_exp_sinh_create_f(void);
QuadratureHandle quad_exp_sinh_create_d(void);
QuadratureHandle quad_exp_sinh_create_l(void);

QuadratureHandle quad_exp_sinh_create_with_params_f(int max_refinements, float tolerance);
QuadratureHandle quad_exp_sinh_create_with_params_d(int max_refinements, double tolerance);
QuadratureHandle quad_exp_sinh_create_with_params_l(int max_refinements, long double tolerance);

// === Shared functions ===
void quad_destroy(QuadratureHandle handle);

// Integration
QuadratureResultF quad_integrate_f(QuadratureHandle handle, IntegrandFunctionF f, void* context);
QuadratureResultD quad_integrate_d(QuadratureHandle handle, IntegrandFunctionD f, void* context);
QuadratureResultL quad_integrate_l(QuadratureHandle handle, IntegrandFunctionL f, void* context);

QuadratureResultF quad_integrate_interval_f(QuadratureHandle handle, IntegrandFunctionF f, void* context, float a, float b);
QuadratureResultD quad_integrate_interval_d(QuadratureHandle handle, IntegrandFunctionD f, void* context, double a, double b);
QuadratureResultL quad_integrate_interval_l(QuadratureHandle handle, IntegrandFunctionL f, void* context, long double a, long double b);

// Information
QuadratureType quad_get_type(QuadratureHandle handle);
QuadraturePrecision quad_get_precision(QuadratureHandle handle);
int quad_get_points(QuadratureHandle handle);

// Abscissa and weights (only for fixed quadratures)
int quad_get_abscissa_weights_f(QuadratureHandle handle, float* abscissa, float* weights, int buffer_size);
int quad_get_abscissa_weights_d(QuadratureHandle handle, double* abscissa, double* weights, int buffer_size);
int quad_get_abscissa_weights_l(QuadratureHandle handle, long double* abscissa, long double* weights, int buffer_size);

// Helper functions
const char* quad_type_to_string(QuadratureType type);
int quad_is_adaptive(QuadratureType type);
int quad_supports_infinite_bounds(QuadratureType type);

#ifdef __cplusplus
}
#endif

#endif // BOOST_QUADRATURE_WRAPPER_HPP
