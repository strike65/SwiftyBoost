//
//  Created by VT on 24.10.25.
//  Copyright © 2025 Volker Thieme. All rights reserved.
//
//  Permission is hereby granted, free of charge, to any person obtaining a copy
//  of this software and associated documentation files (the "Software"), to deal
//  in the Software without restriction, including without limitation the rights
//  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//  copies of the Software, and to permit persons to whom the Software is
//  furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included in
//  all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
//  THE SOFTWARE.
//

#include "../include/Quadrature/bs_quadrature.h"

#include <boost/math/quadrature/exp_sinh.hpp>
#include <boost/math/quadrature/gauss.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/quadrature/sinh_sinh.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>

#include <cmath>
#include <limits>
#include <memory>
#include <type_traits>

// Since Boost 1.78 additional quadrature rules are available.
#if BOOST_VERSION >= 107800
#include <boost/math/quadrature/gauss_hermite.hpp>
#include <boost/math/quadrature/gauss_jacobi.hpp>
#include <boost/math/quadrature/gauss_laguerre.hpp>
#endif

using namespace boost::math::quadrature;

namespace bs::quadrature::detail {

enum class PrecisionTag {
    Float,
    Double,
    LongDouble
};

template <typename Real>
struct Traits;

template <>
struct Traits<float> {
    static constexpr PrecisionTag precision = PrecisionTag::Float;
    using Integrand = IntegrandFunctionF;
    using CResult = QuadratureResultF;
};

template <>
struct Traits<double> {
    static constexpr PrecisionTag precision = PrecisionTag::Double;
    using Integrand = IntegrandFunctionD;
    using CResult = QuadratureResultD;
};

template <>
struct Traits<long double> {
    static constexpr PrecisionTag precision = PrecisionTag::LongDouble;
    using Integrand = IntegrandFunctionL;
    using CResult = QuadratureResultL;
};

class QuadratureAny {
public:
    explicit QuadratureAny(PrecisionTag p) : precision_(p) {}
    virtual ~QuadratureAny() = default;

    virtual QuadratureType get_type() const = 0;
    virtual int get_points() const = 0;

    PrecisionTag precision() const { return precision_; }

private:
    PrecisionTag precision_;
};

template <typename Real>
class QuadratureBase : public QuadratureAny {
public:
    using Integrand = Real (*)(Real, void*);

    struct Result {
        Real result;
        Real error;
        Real l1_norm;
        int iterations;
        int function_calls;
        int converged;
    };

    QuadratureBase() : QuadratureAny(Traits<Real>::precision) {}
    ~QuadratureBase() override = default;

    virtual Result integrate(Integrand f, void* context) = 0;
    virtual Result integrate_interval(Integrand f, void* context, Real a, Real b) = 0;

    virtual bool get_abscissa_weights(Real*, Real*, int) { return false; }
};

template <typename Real>
[[nodiscard]] inline typename QuadratureBase<Real>::Result make_error_result() {
    using Result = typename QuadratureBase<Real>::Result;
    Result r{};
    r.result = Real(0);
    r.error = std::numeric_limits<Real>::infinity();
    r.l1_norm = Real(0);
    r.iterations = 0;
    r.function_calls = 0;
    r.converged = 0;
    return r;
}

template <typename Real>
[[nodiscard]] inline PrecisionTag precision_of() {
    return Traits<Real>::precision;
}

template <typename Real>
[[nodiscard]] inline typename Traits<Real>::CResult to_c_result(const typename QuadratureBase<Real>::Result& src) {
    using CResult = typename Traits<Real>::CResult;
    return CResult{
        static_cast<Real>(src.result),
        static_cast<Real>(src.error),
        static_cast<Real>(src.l1_norm),
        src.iterations,
        src.function_calls,
        src.converged
    };
}

template <typename Real, unsigned Points>
class GaussLegendreImpl final : public QuadratureBase<Real> {
    using Base = QuadratureBase<Real>;
    using Result = typename Base::Result;
    using Integrand = typename Base::Integrand;
    using GaussType = gauss<Real, Points>;

public:
    Result integrate(Integrand f, void* context) override {
        Result result{};
        result.converged = 1;
        result.iterations = 1;
        result.function_calls = static_cast<int>(Points);
        auto lambda = [f, context](Real x) -> Real {
            return f ? f(x, context) : Real(0);
        };
        result.result = GaussType::integrate(lambda);
        result.error = Real(0);
        result.l1_norm = Real(0);
        return result;
    }

    Result integrate_interval(Integrand f, void* context, Real a, Real b) override {
        Result result{};
        result.converged = 1;
        result.iterations = 1;
        result.function_calls = static_cast<int>(Points);
        auto lambda = [f, context](Real x) -> Real {
            return f ? f(x, context) : Real(0);
        };
        result.result = GaussType::integrate(lambda, a, b);
        result.error = Real(0);
        result.l1_norm = Real(0);
        return result;
    }

    QuadratureType get_type() const override { return QUAD_GAUSS_LEGENDRE; }
    int get_points() const override { return static_cast<int>(Points); }

    bool get_abscissa_weights(Real* abscissa, Real* weights, int buffer_size) override {
        if (buffer_size < static_cast<int>(Points)) {
            return false;
        }
        const auto& abs = GaussType::abscissa();
        const auto& wgt = GaussType::weights();
        for (unsigned i = 0; i < Points; ++i) {
            abscissa[i] = abs[i];
            weights[i] = wgt[i];
        }
        return true;
    }
};

#if BOOST_VERSION >= 107800
template <typename Real, unsigned Points>
class GaussHermiteImpl final : public QuadratureBase<Real> {
    using Base = QuadratureBase<Real>;
    using Result = typename Base::Result;
    using Integrand = typename Base::Integrand;
    using GaussType = gauss_hermite<Real, Points>;

public:
    Result integrate(Integrand f, void* context) override {
        Result result{};
        result.converged = 1;
        result.iterations = 1;
        result.function_calls = static_cast<int>(Points);
        auto lambda = [f, context](Real x) -> Real {
            return f ? f(x, context) : Real(0);
        };
        result.result = GaussType::integrate(lambda);
        result.error = Real(0);
        result.l1_norm = Real(0);
        return result;
    }

    Result integrate_interval(Integrand f, void* context, Real, Real) override {
        // Hermite integrates over (-∞, ∞); ignore bounds.
        return integrate(f, context);
    }

    QuadratureType get_type() const override { return QUAD_GAUSS_HERMITE; }
    int get_points() const override { return static_cast<int>(Points); }

    bool get_abscissa_weights(Real* abscissa, Real* weights, int buffer_size) override {
        if (buffer_size < static_cast<int>(Points)) {
            return false;
        }
        const auto& abs = GaussType::abscissa();
        const auto& wgt = GaussType::weights();
        for (unsigned i = 0; i < Points; ++i) {
            abscissa[i] = abs[i];
            weights[i] = wgt[i];
        }
        return true;
    }
};
#endif // BOOST_VERSION >= 107800

#if BOOST_VERSION >= 107800
template <typename Real, unsigned Points>
class GaussLaguerreImpl final : public QuadratureBase<Real> {
    using Base = QuadratureBase<Real>;
    using Result = typename Base::Result;
    using Integrand = typename Base::Integrand;
    using GaussType = gauss_laguerre<Real, Points>;

public:
    explicit GaussLaguerreImpl(Real alpha = Real(0))
        : alpha_(alpha) {}

    Result integrate(Integrand f, void* context) override {
        Result result{};
        result.converged = 1;
        result.iterations = 1;
        result.function_calls = static_cast<int>(Points);
        auto lambda = [f, context](Real x) -> Real {
            return f ? f(x, context) : Real(0);
        };
        if (alpha_ == Real(0)) {
            result.result = GaussType::integrate(lambda);
        } else {
            gauss_laguerre<Real, Points> quad(alpha_);
            result.result = quad.integrate(lambda);
        }
        result.error = Real(0);
        result.l1_norm = Real(0);
        return result;
    }

    Result integrate_interval(Integrand f, void* context, Real, Real) override {
        // Laguerre integrates over [0, ∞); ignore bounds.
        return integrate(f, context);
    }

    QuadratureType get_type() const override { return QUAD_GAUSS_LAGUERRE; }
    int get_points() const override { return static_cast<int>(Points); }

private:
    Real alpha_;
};

template <typename Real, unsigned Points>
class GaussJacobiImpl final : public QuadratureBase<Real> {
    using Base = QuadratureBase<Real>;
    using Result = typename Base::Result;
    using Integrand = typename Base::Integrand;

public:
    GaussJacobiImpl(Real alpha, Real beta)
        : alpha_(alpha), beta_(beta) {}

    Result integrate(Integrand f, void* context) override {
        return integrate_interval(f, context, Real(-1), Real(1));
    }

    Result integrate_interval(Integrand f, void* context, Real a, Real b) override {
        Result result{};
        result.converged = 1;
        result.iterations = 1;
        result.function_calls = static_cast<int>(Points);
        gauss_jacobi<Real, Points> quad(alpha_, beta_);
        auto lambda = [f, context](Real x) -> Real {
            return f ? f(x, context) : Real(0);
        };
        result.result = quad.integrate(lambda, a, b);
        result.error = Real(0);
        result.l1_norm = Real(0);
        return result;
    }

    QuadratureType get_type() const override { return QUAD_GAUSS_JACOBI; }
    int get_points() const override { return static_cast<int>(Points); }

private:
    Real alpha_;
    Real beta_;
};
#endif // BOOST_VERSION >= 107800

template <typename Real, unsigned Points>
class GaussKronrodImpl final : public QuadratureBase<Real> {
    using Base = QuadratureBase<Real>;
    using Result = typename Base::Result;
    using Integrand = typename Base::Integrand;
    using GaussKronrodType = gauss_kronrod<Real, Points>;

public:
    Result integrate(Integrand f, void* context) override {
        return integrate_interval(f, context, Real(-1), Real(1));
    }

    Result integrate_interval(Integrand f, void* context, Real a, Real b) override {
        Result result{};
        result.converged = 1;
        result.iterations = 1;
        result.function_calls = static_cast<int>(Points);
        auto lambda = [f, context](Real x) -> Real {
            return f ? f(x, context) : Real(0);
        };
        result.result = GaussKronrodType::integrate(lambda, a, b);
        result.error = Real(0);
        result.l1_norm = Real(0);
        return result;
    }

    QuadratureType get_type() const override { return QUAD_GAUSS_KRONROD; }
    int get_points() const override { return static_cast<int>(Points); }
};

template <typename Real>
class TanhSinhImpl final : public QuadratureBase<Real> {
    using Base = QuadratureBase<Real>;
    using Result = typename Base::Result;
    using Integrand = typename Base::Integrand;

public:
    TanhSinhImpl(int max_ref = 10, Real tol = Real(1e-9))
        : integrator_(max_ref),
          max_refinements_(max_ref),
          tolerance_(tol) {}

    Result integrate(Integrand f, void* context) override {
        return integrate_interval(f, context, Real(-1), Real(1));
    }

    Result integrate_interval(Integrand f, void* context, Real a, Real b) override {
        Result result{};
        result.converged = 1;
        result.iterations = 0;
        result.function_calls = 0;
        auto lambda = [f, context, &result](Real x) -> Real {
            result.function_calls++;
            return f ? f(x, context) : Real(0);
        };
        result.result = integrator_.integrate(lambda, a, b, tolerance_);
        result.error = Real(0);
        result.l1_norm = Real(0);
        return result;
    }

    QuadratureType get_type() const override { return QUAD_TANH_SINH; }
    int get_points() const override { return -1; }

private:
    tanh_sinh<Real> integrator_;
    int max_refinements_;
    Real tolerance_;
};

template <typename Real>
class SinhSinhImpl final : public QuadratureBase<Real> {
    using Base = QuadratureBase<Real>;
    using Result = typename Base::Result;
    using Integrand = typename Base::Integrand;

public:
    SinhSinhImpl(int max_ref = 10, Real tol = Real(1e-9))
        : integrator_(max_ref),
          max_refinements_(max_ref),
          tolerance_(tol) {}

    Result integrate(Integrand f, void* context) override {
        Result result{};
        result.converged = 1;
        result.iterations = 0;
        result.function_calls = 0;
        auto lambda = [f, context, &result](Real x) -> Real {
            result.function_calls++;
            return f ? f(x, context) : Real(0);
        };
        result.result = integrator_.integrate(lambda, tolerance_);
        result.error = Real(0);
        result.l1_norm = Real(0);
        return result;
    }

    Result integrate_interval(Integrand f, void* context, Real, Real) override {
        // Sinh-Sinh integrates over (-∞, ∞); ignore bounds.
        return integrate(f, context);
    }

    QuadratureType get_type() const override { return QUAD_SINH_SINH; }
    int get_points() const override { return -1; }

private:
    sinh_sinh<Real> integrator_;
    int max_refinements_;
    Real tolerance_;
};

template <typename Real>
class ExpSinhImpl final : public QuadratureBase<Real> {
    using Base = QuadratureBase<Real>;
    using Result = typename Base::Result;
    using Integrand = typename Base::Integrand;

public:
    ExpSinhImpl(int max_ref = 10, Real tol = Real(1e-9))
        : integrator_(max_ref),
          max_refinements_(max_ref),
          tolerance_(tol) {}

    Result integrate(Integrand f, void* context) override {
        Result result{};
        result.converged = 1;
        result.iterations = 0;
        result.function_calls = 0;
        auto lambda = [f, context, &result](Real x) -> Real {
            result.function_calls++;
            return f ? f(x, context) : Real(0);
        };
        result.result = integrator_.integrate(lambda, tolerance_);
        result.error = Real(0);
        result.l1_norm = Real(0);
        return result;
    }

    Result integrate_interval(Integrand f, void* context, Real a, Real) override {
        if (a != Real(0)) {
            auto shifted = [f, context, a](Real x) -> Real {
                return f ? f(x + a, context) : Real(0);
            };
            Result result{};
            result.converged = 1;
            result.iterations = 0;
            result.function_calls = 0;
            auto lambda = [shifted, &result](Real x) -> Real {
                result.function_calls++;
                return shifted(x);
            };
            result.result = integrator_.integrate(lambda, tolerance_);
            result.error = Real(0);
            result.l1_norm = Real(0);
            return result;
        }
        return integrate(f, context);
    }

    QuadratureType get_type() const override { return QUAD_EXP_SINH; }
    int get_points() const override { return -1; }

private:
    exp_sinh<Real> integrator_;
    int max_refinements_;
    Real tolerance_;
};

template <typename Real>
QuadratureBase<Real>* create_gauss(int points) {
    switch (points) {
        case 7:  return new GaussLegendreImpl<Real, 7>();
        case 10: return new GaussLegendreImpl<Real, 10>();
        case 15: return new GaussLegendreImpl<Real, 15>();
        case 20: return new GaussLegendreImpl<Real, 20>();
        case 25: return new GaussLegendreImpl<Real, 25>();
        case 30: return new GaussLegendreImpl<Real, 30>();
        case 40: return new GaussLegendreImpl<Real, 40>();
        case 50: return new GaussLegendreImpl<Real, 50>();
        case 60: return new GaussLegendreImpl<Real, 60>();
        case 70: return new GaussLegendreImpl<Real, 70>();
        case 80: return new GaussLegendreImpl<Real, 80>();
        case 90: return new GaussLegendreImpl<Real, 90>();
        case 100: return new GaussLegendreImpl<Real, 100>();
        default: return nullptr;
    }
}

#if BOOST_VERSION >= 107800
template <typename Real>
QuadratureBase<Real>* create_gauss_hermite(int points) {
    switch (points) {
        case 10: return new GaussHermiteImpl<Real, 10>();
        case 15: return new GaussHermiteImpl<Real, 15>();
        case 20: return new GaussHermiteImpl<Real, 20>();
        case 25: return new GaussHermiteImpl<Real, 25>();
        case 30: return new GaussHermiteImpl<Real, 30>();
        case 40: return new GaussHermiteImpl<Real, 40>();
        case 50: return new GaussHermiteImpl<Real, 50>();
        default: return nullptr;
    }
}

template <typename Real>
QuadratureBase<Real>* create_gauss_laguerre(int points, Real alpha, bool use_alpha) {
    switch (points) {
        case 10: return use_alpha ? static_cast<QuadratureBase<Real>*>(new GaussLaguerreImpl<Real, 10>(alpha)) : new GaussLaguerreImpl<Real, 10>();
        case 15: return use_alpha ? static_cast<QuadratureBase<Real>*>(new GaussLaguerreImpl<Real, 15>(alpha)) : new GaussLaguerreImpl<Real, 15>();
        case 20: return use_alpha ? static_cast<QuadratureBase<Real>*>(new GaussLaguerreImpl<Real, 20>(alpha)) : new GaussLaguerreImpl<Real, 20>();
        case 25: return use_alpha ? static_cast<QuadratureBase<Real>*>(new GaussLaguerreImpl<Real, 25>(alpha)) : new GaussLaguerreImpl<Real, 25>();
        case 30: return use_alpha ? static_cast<QuadratureBase<Real>*>(new GaussLaguerreImpl<Real, 30>(alpha)) : new GaussLaguerreImpl<Real, 30>();
        case 40: return use_alpha ? static_cast<QuadratureBase<Real>*>(new GaussLaguerreImpl<Real, 40>(alpha)) : new GaussLaguerreImpl<Real, 40>();
        case 50: return use_alpha ? static_cast<QuadratureBase<Real>*>(new GaussLaguerreImpl<Real, 50>(alpha)) : new GaussLaguerreImpl<Real, 50>();
        default: return nullptr;
    }
}

template <typename Real>
QuadratureBase<Real>* create_gauss_jacobi(int points, Real alpha, Real beta) {
    switch (points) {
        case 10: return new GaussJacobiImpl<Real, 10>(alpha, beta);
        case 15: return new GaussJacobiImpl<Real, 15>(alpha, beta);
        case 20: return new GaussJacobiImpl<Real, 20>(alpha, beta);
        case 30: return new GaussJacobiImpl<Real, 30>(alpha, beta);
        case 50: return new GaussJacobiImpl<Real, 50>(alpha, beta);
        default: return nullptr;
    }
}
#endif // BOOST_VERSION >= 107800

template <typename Real>
QuadratureBase<Real>* create_gauss_kronrod(int points) {
    switch (points) {
        case 15: return new GaussKronrodImpl<Real, 15>();
        case 21: return new GaussKronrodImpl<Real, 21>();
        case 31: return new GaussKronrodImpl<Real, 31>();
        case 41: return new GaussKronrodImpl<Real, 41>();
        case 51: return new GaussKronrodImpl<Real, 51>();
        case 61: return new GaussKronrodImpl<Real, 61>();
        default: return nullptr;
    }
}

template <typename Real>
QuadratureBase<Real>* create_tanh_sinh(int max_refinements, Real tolerance, bool custom) {
    if (custom) {
        return new TanhSinhImpl<Real>(max_refinements, tolerance);
    }
    return new TanhSinhImpl<Real>();
}

template <typename Real>
QuadratureBase<Real>* create_sinh_sinh(int max_refinements, Real tolerance, bool custom) {
    if (custom) {
        return new SinhSinhImpl<Real>(max_refinements, tolerance);
    }
    return new SinhSinhImpl<Real>();
}

template <typename Real>
QuadratureBase<Real>* create_exp_sinh(int max_refinements, Real tolerance, bool custom) {
    if (custom) {
        return new ExpSinhImpl<Real>(max_refinements, tolerance);
    }
    return new ExpSinhImpl<Real>();
}

template <typename Real>
QuadratureBase<Real>* as_impl(QuadratureHandle handle) {
    return dynamic_cast<QuadratureBase<Real>*>(static_cast<QuadratureAny*>(handle));
}

inline QuadratureAny* as_any(QuadratureHandle handle) {
    return static_cast<QuadratureAny*>(handle);
}

template <typename Real>
typename QuadratureBase<Real>::Result integrate_impl(QuadratureHandle handle,
                                                     typename QuadratureBase<Real>::Integrand f,
                                                     void* context) {
    auto* impl = as_impl<Real>(handle);
    if (!impl) {
        return make_error_result<Real>();
    }
    return impl->integrate(f, context);
}

template <typename Real>
typename QuadratureBase<Real>::Result integrate_interval_impl(QuadratureHandle handle,
                                                              typename QuadratureBase<Real>::Integrand f,
                                                              void* context,
                                                              Real a,
                                                              Real b) {
    auto* impl = as_impl<Real>(handle);
    if (!impl) {
        return make_error_result<Real>();
    }
    return impl->integrate_interval(f, context, a, b);
}

template <typename Real>
int get_abscissa_impl(QuadratureHandle handle, Real* abscissa, Real* weights, int buffer_size) {
    auto* impl = as_impl<Real>(handle);
    if (!impl) {
        return 0;
    }
    return impl->get_abscissa_weights(abscissa, weights, buffer_size) ? 1 : 0;
}

inline QuadraturePrecision to_c_precision(PrecisionTag tag) {
    switch (tag) {
        case PrecisionTag::Float: return QUAD_PRECISION_FLOAT;
        case PrecisionTag::Double: return QUAD_PRECISION_DOUBLE;
        case PrecisionTag::LongDouble: return QUAD_PRECISION_LONG_DOUBLE;
    }
    return QUAD_PRECISION_DOUBLE;
}

template <typename Real, typename Factory>
QuadratureHandle create_handle(Factory&& factory) {
    std::unique_ptr<QuadratureBase<Real>> ptr(factory());
    if (!ptr) {
        return nullptr;
    }
    return ptr.release();
}

} // namespace bs::quadrature::detail

using namespace bs::quadrature::detail;

// === Factory implementations =================================================

QuadratureHandle quad_gauss_create_f(int points) {
    return create_handle<float>([&]() { return create_gauss<float>(points); });
}

QuadratureHandle quad_gauss_create_d(int points) {
    return create_handle<double>([&]() { return create_gauss<double>(points); });
}

QuadratureHandle quad_gauss_create_l(int points) {
    return create_handle<long double>([&]() { return create_gauss<long double>(points); });
}

#if BOOST_VERSION >= 107800
QuadratureHandle quad_gauss_hermite_create_f(int points) {
    return create_handle<float>([&]() { return create_gauss_hermite<float>(points); });
}

QuadratureHandle quad_gauss_hermite_create_d(int points) {
    return create_handle<double>([&]() { return create_gauss_hermite<double>(points); });
}

QuadratureHandle quad_gauss_hermite_create_l(int points) {
    return create_handle<long double>([&]() { return create_gauss_hermite<long double>(points); });
}
#else
QuadratureHandle quad_gauss_hermite_create_f(int) { return nullptr; }
QuadratureHandle quad_gauss_hermite_create_d(int) { return nullptr; }
QuadratureHandle quad_gauss_hermite_create_l(int) { return nullptr; }
#endif

#if BOOST_VERSION >= 107800
QuadratureHandle quad_gauss_laguerre_create_f(int points) {
    return create_handle<float>([&]() { return create_gauss_laguerre<float>(points, float(0), false); });
}

QuadratureHandle quad_gauss_laguerre_create_d(int points) {
    return create_handle<double>([&]() { return create_gauss_laguerre<double>(points, 0.0, false); });
}

QuadratureHandle quad_gauss_laguerre_create_l(int points) {
    return create_handle<long double>([&]() { return create_gauss_laguerre<long double>(points, 0.0L, false); });
}

QuadratureHandle quad_gauss_laguerre_create_alpha_f(int points, float alpha) {
    return create_handle<float>([&]() { return create_gauss_laguerre<float>(points, alpha, true); });
}

QuadratureHandle quad_gauss_laguerre_create_alpha_d(int points, double alpha) {
    return create_handle<double>([&]() { return create_gauss_laguerre<double>(points, alpha, true); });
}

QuadratureHandle quad_gauss_laguerre_create_alpha_l(int points, long double alpha) {
    return create_handle<long double>([&]() { return create_gauss_laguerre<long double>(points, alpha, true); });
}
#else
QuadratureHandle quad_gauss_laguerre_create_f(int) { return nullptr; }
QuadratureHandle quad_gauss_laguerre_create_d(int) { return nullptr; }
QuadratureHandle quad_gauss_laguerre_create_l(int) { return nullptr; }
QuadratureHandle quad_gauss_laguerre_create_alpha_f(int, float) { return nullptr; }
QuadratureHandle quad_gauss_laguerre_create_alpha_d(int, double) { return nullptr; }
QuadratureHandle quad_gauss_laguerre_create_alpha_l(int, long double) { return nullptr; }
#endif

#if BOOST_VERSION >= 107800
QuadratureHandle quad_gauss_jacobi_create_f(int points, float alpha, float beta) {
    return create_handle<float>([&]() { return create_gauss_jacobi<float>(points, alpha, beta); });
}

QuadratureHandle quad_gauss_jacobi_create_d(int points, double alpha, double beta) {
    return create_handle<double>([&]() { return create_gauss_jacobi<double>(points, alpha, beta); });
}

QuadratureHandle quad_gauss_jacobi_create_l(int points, long double alpha, long double beta) {
    return create_handle<long double>([&]() { return create_gauss_jacobi<long double>(points, alpha, beta); });
}
#else
QuadratureHandle quad_gauss_jacobi_create_f(int, float, float) { return nullptr; }
QuadratureHandle quad_gauss_jacobi_create_d(int, double, double) { return nullptr; }
QuadratureHandle quad_gauss_jacobi_create_l(int, long double, long double) { return nullptr; }
#endif

QuadratureHandle quad_gauss_kronrod_create_f(int points) {
    return create_handle<float>([&]() { return create_gauss_kronrod<float>(points); });
}

QuadratureHandle quad_gauss_kronrod_create_d(int points) {
    return create_handle<double>([&]() { return create_gauss_kronrod<double>(points); });
}

QuadratureHandle quad_gauss_kronrod_create_l(int points) {
    return create_handle<long double>([&]() { return create_gauss_kronrod<long double>(points); });
}

QuadratureHandle quad_tanh_sinh_create_f(void) {
    return create_handle<float>([&]() { return create_tanh_sinh<float>(10, float(1e-9), false); });
}

QuadratureHandle quad_tanh_sinh_create_d(void) {
    return create_handle<double>([&]() { return create_tanh_sinh<double>(10, 1e-9, false); });
}

QuadratureHandle quad_tanh_sinh_create_l(void) {
    return create_handle<long double>([&]() { return create_tanh_sinh<long double>(10, 1e-9L, false); });
}

QuadratureHandle quad_tanh_sinh_create_with_params_f(int max_ref, float tol) {
    return create_handle<float>([&]() { return create_tanh_sinh<float>(max_ref, tol, true); });
}

QuadratureHandle quad_tanh_sinh_create_with_params_d(int max_ref, double tol) {
    return create_handle<double>([&]() { return create_tanh_sinh<double>(max_ref, tol, true); });
}

QuadratureHandle quad_tanh_sinh_create_with_params_l(int max_ref, long double tol) {
    return create_handle<long double>([&]() { return create_tanh_sinh<long double>(max_ref, tol, true); });
}

QuadratureHandle quad_sinh_sinh_create_f(void) {
    return create_handle<float>([&]() { return create_sinh_sinh<float>(10, float(1e-9), false); });
}

QuadratureHandle quad_sinh_sinh_create_d(void) {
    return create_handle<double>([&]() { return create_sinh_sinh<double>(10, 1e-9, false); });
}

QuadratureHandle quad_sinh_sinh_create_l(void) {
    return create_handle<long double>([&]() { return create_sinh_sinh<long double>(10, 1e-9L, false); });
}

QuadratureHandle quad_sinh_sinh_create_with_params_f(int max_ref, float tol) {
    return create_handle<float>([&]() { return create_sinh_sinh<float>(max_ref, tol, true); });
}

QuadratureHandle quad_sinh_sinh_create_with_params_d(int max_ref, double tol) {
    return create_handle<double>([&]() { return create_sinh_sinh<double>(max_ref, tol, true); });
}

QuadratureHandle quad_sinh_sinh_create_with_params_l(int max_ref, long double tol) {
    return create_handle<long double>([&]() { return create_sinh_sinh<long double>(max_ref, tol, true); });
}

QuadratureHandle quad_exp_sinh_create_f(void) {
    return create_handle<float>([&]() { return create_exp_sinh<float>(10, float(1e-9), false); });
}

QuadratureHandle quad_exp_sinh_create_d(void) {
    return create_handle<double>([&]() { return create_exp_sinh<double>(10, 1e-9, false); });
}

QuadratureHandle quad_exp_sinh_create_l(void) {
    return create_handle<long double>([&]() { return create_exp_sinh<long double>(10, 1e-9L, false); });
}

QuadratureHandle quad_exp_sinh_create_with_params_f(int max_ref, float tol) {
    return create_handle<float>([&]() { return create_exp_sinh<float>(max_ref, tol, true); });
}

QuadratureHandle quad_exp_sinh_create_with_params_d(int max_ref, double tol) {
    return create_handle<double>([&]() { return create_exp_sinh<double>(max_ref, tol, true); });
}

QuadratureHandle quad_exp_sinh_create_with_params_l(int max_ref, long double tol) {
    return create_handle<long double>([&]() { return create_exp_sinh<long double>(max_ref, tol, true); });
}

// === Lifetime ================================================================

void quad_destroy(QuadratureHandle handle) {
    delete as_any(handle);
}

// === Integration =============================================================

QuadratureResultF quad_integrate_f(QuadratureHandle handle, IntegrandFunctionF f, void* context) {
    auto raw = integrate_impl<float>(handle, f, context);
    return to_c_result<float>(raw);
}

QuadratureResultD quad_integrate_d(QuadratureHandle handle, IntegrandFunctionD f, void* context) {
    auto raw = integrate_impl<double>(handle, f, context);
    return to_c_result<double>(raw);
}

QuadratureResultL quad_integrate_l(QuadratureHandle handle, IntegrandFunctionL f, void* context) {
    auto raw = integrate_impl<long double>(handle, f, context);
    return to_c_result<long double>(raw);
}

QuadratureResultF quad_integrate_interval_f(QuadratureHandle handle,
                                            IntegrandFunctionF f,
                                            void* context,
                                            float a,
                                            float b) {
    auto raw = integrate_interval_impl<float>(handle, f, context, a, b);
    return to_c_result<float>(raw);
}

QuadratureResultD quad_integrate_interval_d(QuadratureHandle handle,
                                            IntegrandFunctionD f,
                                            void* context,
                                            double a,
                                            double b) {
    auto raw = integrate_interval_impl<double>(handle, f, context, a, b);
    return to_c_result<double>(raw);
}

QuadratureResultL quad_integrate_interval_l(QuadratureHandle handle,
                                            IntegrandFunctionL f,
                                            void* context,
                                            long double a,
                                            long double b) {
    auto raw = integrate_interval_impl<long double>(handle, f, context, a, b);
    return to_c_result<long double>(raw);
}

// === Metadata ================================================================

QuadratureType quad_get_type(QuadratureHandle handle) {
    auto* any = as_any(handle);
    return any ? any->get_type() : QUAD_GAUSS_LEGENDRE;
}

QuadraturePrecision quad_get_precision(QuadratureHandle handle) {
    auto* any = as_any(handle);
    return any ? to_c_precision(any->precision()) : QUAD_PRECISION_DOUBLE;
}

int quad_get_points(QuadratureHandle handle) {
    auto* any = as_any(handle);
    return any ? any->get_points() : 0;
}

int quad_get_abscissa_weights_f(QuadratureHandle handle, float* abscissa, float* weights, int buffer_size) {
    return get_abscissa_impl<float>(handle, abscissa, weights, buffer_size);
}

int quad_get_abscissa_weights_d(QuadratureHandle handle, double* abscissa, double* weights, int buffer_size) {
    return get_abscissa_impl<double>(handle, abscissa, weights, buffer_size);
}

int quad_get_abscissa_weights_l(QuadratureHandle handle, long double* abscissa, long double* weights, int buffer_size) {
    return get_abscissa_impl<long double>(handle, abscissa, weights, buffer_size);
}

// === Helper utilities ========================================================

const char* quad_type_to_string(QuadratureType type) {
    switch (type) {
        case QUAD_GAUSS_LEGENDRE: return "gauss_legendre";
        case QUAD_GAUSS_HERMITE: return "gauss_hermite";
        case QUAD_GAUSS_LAGUERRE: return "gauss_laguerre";
        case QUAD_GAUSS_JACOBI: return "gauss_jacobi";
        case QUAD_GAUSS_KRONROD: return "gauss_kronrod";
        case QUAD_TANH_SINH: return "tanh_sinh";
        case QUAD_SINH_SINH: return "sinh_sinh";
        case QUAD_EXP_SINH: return "exp_sinh";
    }
    return "unknown";
}

int quad_is_adaptive(QuadratureType type) {
    switch (type) {
        case QUAD_TANH_SINH:
        case QUAD_SINH_SINH:
        case QUAD_EXP_SINH:
            return 1;
        default:
            return 0;
    }
}

int quad_supports_infinite_bounds(QuadratureType type) {
    switch (type) {
        case QUAD_GAUSS_HERMITE:
        case QUAD_SINH_SINH:
        case QUAD_EXP_SINH:
            return 1;
        default:
            return 0;
    }
}
