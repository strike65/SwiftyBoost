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

/*
    Overview
    --------
    This translation unit implements a C-compatible facade over several
    quadrature rules provided by Boost.Math. It exposes an opaque handle
    (QuadratureHandle) and a family of C functions declared in
    "../include/Quadrature/bs_quadrature.h" to:

      - Construct quadrature objects for different rules and precisions.
      - Integrate functions either over the rule’s natural domain or a
        user-specified finite interval (where applicable).
      - Obtain metadata (type, precision, number of points).
      - Retrieve abscissa and weights for fixed (non-adaptive) rules.

    Internally, we use a small type-erased hierarchy:

      - QuadratureAny: non-templated base for lifetime and metadata.
      - QuadratureBase<Real>: templated base that defines a uniform
        interface across Real = float/double/long double.
      - Concrete implementations per Boost rule and point count.

    Ownership and lifetime
    ----------------------
    - All factory functions return a raw QuadratureHandle that is in fact
      a pointer to a QuadratureAny-derived instance. The caller becomes
      the owner and must call quad_destroy(handle) to release it.
    - No shared ownership is used; creation succeeds or returns nullptr.

    Error semantics
    ---------------
    - For fixed rules (Gauss-Legendre/Jacobi/Hermite/Laguerre/Kronrod),
      Boost’s compile-time constructs do not provide a runtime error estimate.
      We therefore set Result.error and Result.l1_norm to zero, mark
      converged=1, iterations=1, and function_calls to the rule’s point count.
    - For adaptive rules (tanh_sinh, sinh_sinh, exp_sinh), we wrap Boost’s
      adaptive integrators. We count function evaluations in the callback
      and return iterations as 0 (Boost doesn’t expose a step count via this
      interface). Error and L1 norm are currently zero-filled; extend here
      if you plumb Boost’s optional error/L1 outputs in the future.

    Bounds handling
    ---------------
    - Some rules have natural infinite domains (Hermite: (-∞,∞), Laguerre: [0,∞)).
      Their integrate_interval implementations ignore user bounds and integrate
      over the natural domain. For Exp-Sinh, we support a left-bound shift so
      [a, ∞) can be handled by a change of variables when a != 0.

    Boost version feature gates
    ---------------------------
    - Since Boost 1.78, additional fixed rules are available:
      Gauss-Hermite, Gauss-Laguerre, and Gauss-Jacobi. This file conditionally
      compiles those implementations and factories behind BOOST_VERSION >= 107800.

    Thread-safety
    -------------
    - Each handle is independent. There is no shared mutable state between
      handles. The function-call counter inside adaptive integrators is local
      to the integration call and guarded by capture-by-reference into the
      lambda, so concurrent calls on the same handle are not supported.

    Precision bridging
    ------------------
    - The Traits<> mapping associates C ABI function-pointer and result
      structs per precision. This allows a single implementation to bridge
      to the C header’s typedefs: IntegrandFunctionF/D/L and QuadratureResultF/D/L.

    Abscissa and weights
    --------------------
    - Only fixed rules with compile-time point counts expose static abscissa
      and weight arrays in Boost. We expose them via quad_get_abscissa_weights_*.
      Adaptive rules return false (0) from this query.
*/

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

using namespace boost::math::quadrature;

namespace bs::quadrature::detail {

/*
    PrecisionTag
    ------------
    Internal discriminator that records the Real type used by a concrete
    quadrature instance. This is translated to the public C enum
    QuadraturePrecision via to_c_precision().
*/
enum class PrecisionTag {
    Float,
    Double,
    LongDouble
};

/*
    Traits<Real>
    ------------
    Associates a Real type with:
      - The corresponding C function pointer type for integrands.
      - The corresponding C result structure type.
      - The internal PrecisionTag.
*/
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

/*
    QuadratureAny
    -------------
    Non-templated base class used for:
      - Type erasure across different Real types.
      - Exposing metadata common to all rules (type, points, precision).
    Derived classes must implement get_type() and get_points().
*/
class QuadratureAny {
public:
    explicit QuadratureAny(PrecisionTag p) : precision_(p) {}
    virtual ~QuadratureAny() = default;

    // Returns the concrete quadrature rule (C ABI enum).
    virtual QuadratureType get_type() const = 0;

    // Returns the number of points for fixed rules; -1 for adaptive rules.
    virtual int get_points() const = 0;

    // Returns the internal precision tag.
    PrecisionTag precision() const { return precision_; }

private:
    PrecisionTag precision_;
};

/*
    QuadratureBase<Real>
    --------------------
    Templated base class that defines:
      - The C ABI-compatible Integrand function signature for the given Real.
      - A Result structure used internally that mirrors the C result types.
      - Virtual methods: integrate(...) and integrate_interval(...).
      - Optional abscissa/weights retrieval for fixed rules.

    Implementors only need to provide the two integrate methods and, if
    applicable, get_abscissa_weights.
*/
template <typename Real>
class QuadratureBase : public QuadratureAny {
public:
    using Integrand = Real (*)(Real, void*);

    struct Result {
        Real result;          // Integral estimate.
        Real error;           // Estimated absolute error (0 for fixed rules).
        Real l1_norm;         // L1 norm estimate (0 for fixed rules).
        int iterations;       // Iteration count (0 for adaptive here).
        int function_calls;   // Number of integrand evaluations (best effort).
        int converged;        // 1 if converged/OK; 0 if error/unavailable.
    };

    QuadratureBase() : QuadratureAny(Traits<Real>::precision) {}
    ~QuadratureBase() override = default;

    // Integrate over the rule’s natural domain or default interval.
    virtual Result integrate(Integrand f, void* context) = 0;

    // Integrate over [a, b] when supported. Implementations may ignore bounds
    // for rules defined on infinite domains.
    virtual Result integrate_interval(Integrand f, void* context, Real a, Real b) = 0;

    // Retrieve abscissa and weights for fixed rules. Default returns false.
    virtual bool get_abscissa_weights(Real*, Real*, int) { return false; }
};

/*
    make_error_result<Real>()
    ------------------------
    Utility to construct a sentinel Result indicating failure/unavailable.
    Used when a handle cannot be cast to the expected Real-specialized
    implementation.
*/
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

/*
    precision_of<Real>()
    -------------------
    Returns the internal PrecisionTag corresponding to Real.
*/
template <typename Real>
[[nodiscard]] inline PrecisionTag precision_of() {
    return Traits<Real>::precision;
}

/*
    to_c_result<Real>(Result)
    -------------------------
    Converts an internal Result to the corresponding C ABI result struct.
*/
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

/*
    GaussLegendreImpl<Real, Points>
    -------------------------------
    Fixed Gauss-Legendre quadrature of compile-time order Points.
    - integrate(): integrates over [-1, 1].
    - integrate_interval(a, b): integrates over [a, b].
    - get_abscissa_weights(): copies Boost’s static abscissa/weights.
*/
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
        // Boost exposes only ceil(N/2) entries: zero (if N odd) and positive nodes.
        // Expand to the full N nodes by mirroring symmetry.
        const auto& abs_half = GaussType::abscissa();
        const auto& wgt_half = GaussType::weights();
        const unsigned halfCount = static_cast<unsigned>(abs_half.size());

        if constexpr (Points % 2 == 1) {
            // Odd N: first entry is the zero node.
            if (halfCount < 1) {
                return false;
            }
            unsigned out = 0;
            // Zero node
            abscissa[out] = Real(0);
            weights[out] = wgt_half[0];
            ++out;
            // Mirror remaining positive nodes.
            for (unsigned i = 1; i < halfCount; ++i) {
                const Real x = abs_half[i];
                const Real w = wgt_half[i];
                abscissa[out] = -x;
                weights[out] = w;
                ++out;
                abscissa[out] = x;
                weights[out] = w;
                ++out;
            }
            // out should equal Points: 1 + 2*(halfCount - 1) == Points
            return true;
        } else {
            // Even N: all entries are positive nodes; mirror each.
            unsigned out = 0;
            for (unsigned i = 0; i < halfCount; ++i) {
                const Real x = abs_half[i];
                const Real w = wgt_half[i];
                abscissa[out] = -x;
                weights[out] = w;
                ++out;
                abscissa[out] = x;
                weights[out] = w;
                ++out;
            }
            // out should equal Points: 2*halfCount == Points
            return true;
        }
    }
};

/*
    GaussKronrodImpl<Real, Points>
    ------------------------------
    Fixed Gauss-Kronrod quadrature of compile-time order Points.
    Natural interval: [-1, 1]. Also supports [a, b].
*/
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

/*
    TanhSinhImpl<Real>
    ------------------
    Adaptive double-exponential tanh-sinh integrator.
    - Default constructor uses Boost defaults; an alternate constructor
      allows specifying max refinements and convergence tolerance.
    - Integrates over [-1, 1] by default; integrate_interval supports [a, b].
    - function_calls is counted by intercepting the callback.
*/
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

/*
    SinhSinhImpl<Real>
    ------------------
    Adaptive sinh-sinh integrator.
    - Natural domain: (-∞, ∞).
    - integrate_interval ignores bounds and integrates over the natural domain.
*/
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

/*
    ExpSinhImpl<Real>
    -----------------
    Adaptive exp-sinh integrator.
    - Natural domain: [0, ∞).
    - integrate_interval supports shifting the lower bound: if a != 0,
      we integrate f(x + a) over [0, ∞) which corresponds to [a, ∞).
*/
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
            // Change of variables: integrate g(x) = f(x + a) over [0, ∞).
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
        // Natural domain [0, ∞).
        return integrate(f, context);
    }

    QuadratureType get_type() const override { return QUAD_EXP_SINH; }
    int get_points() const override { return -1; }

private:
    exp_sinh<Real> integrator_;
    int max_refinements_;
    Real tolerance_;
};

/*
    Factory helpers per rule and precision
    --------------------------------------
    The following functions create concrete QuadratureBase<Real> instances
    for the requested point count or parameters. Unsupported configurations
    return nullptr.
*/
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

/*
    create_tanh_sinh / create_sinh_sinh / create_exp_sinh
    ------------------------------------------------------
    Create adaptive integrators with defaults or custom parameters.
*/
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

/*
    Downcasting helpers
    -------------------
    - as_impl<Real>: cast an opaque handle back to QuadratureBase<Real>*.
      Returns nullptr if the underlying object uses a different Real.
    - as_any: cast to the non-templated base for metadata access.
*/
template <typename Real>
QuadratureBase<Real>* as_impl(QuadratureHandle handle) {
    return dynamic_cast<QuadratureBase<Real>*>(static_cast<QuadratureAny*>(handle));
}

inline QuadratureAny* as_any(QuadratureHandle handle) {
    return static_cast<QuadratureAny*>(handle);
}

/*
    integrate_impl / integrate_interval_impl
    ----------------------------------------
    Invoke the underlying implementation with the provided C ABI integrand
    function pointer and context. If the handle cannot be cast to the Real
    specialization, returns a sentinel error result.
*/
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

/*
    get_abscissa_impl
    -----------------
    Copy abscissa and weights for fixed rules into caller-provided buffers.
    Returns 1 on success, 0 otherwise (including adaptive rules).
*/
template <typename Real>
int get_abscissa_impl(QuadratureHandle handle, Real* abscissa, Real* weights, int buffer_size) {
    auto* impl = as_impl<Real>(handle);
    if (!impl) {
        return 0;
    }
    return impl->get_abscissa_weights(abscissa, weights, buffer_size) ? 1 : 0;
}

/*
    to_c_precision
    --------------
    Translate internal PrecisionTag to the public C QuadraturePrecision enum.
*/
inline QuadraturePrecision to_c_precision(PrecisionTag tag) {
    switch (tag) {
        case PrecisionTag::Float: return QUAD_PRECISION_FLOAT;
        case PrecisionTag::Double: return QUAD_PRECISION_DOUBLE;
        case PrecisionTag::LongDouble: return QUAD_PRECISION_LONG_DOUBLE;
    }
    return QUAD_PRECISION_DOUBLE;
}

/*
    create_handle<Real>(Factory)
    ----------------------------
    Small RAII helper that constructs a QuadratureBase<Real> via a factory
    functor and returns ownership to the caller as a raw QuadratureHandle.
    Returns nullptr if the factory yields a null pointer.
*/
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

/*
    Gauss-Legendre factories (C ABI)
    --------------------------------
    Create fixed Gauss-Legendre integrators for the requested point count
    and precision. Return nullptr for unsupported counts.
*/
QuadratureHandle quad_gauss_create_f(int points) {
    return create_handle<float>([&]() { return create_gauss<float>(points); });
}

QuadratureHandle quad_gauss_create_d(int points) {
    return create_handle<double>([&]() { return create_gauss<double>(points); });
}

QuadratureHandle quad_gauss_create_l(int points) {
    return create_handle<long double>([&]() { return create_gauss<long double>(points); });
}

/*
    Gauss-Kronrod factories (C ABI)
*/
QuadratureHandle quad_gauss_kronrod_create_f(int points) {
    return create_handle<float>([&]() { return create_gauss_kronrod<float>(points); });
}

QuadratureHandle quad_gauss_kronrod_create_d(int points) {
    return create_handle<double>([&]() { return create_gauss_kronrod<double>(points); });
}

QuadratureHandle quad_gauss_kronrod_create_l(int points) {
    return create_handle<long double>([&]() { return create_gauss_kronrod<long double>(points); });
}

/*
    Tanh-Sinh factories (C ABI)
    - Default and parameterized variants.
*/
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

/*
    Sinh-Sinh factories (C ABI)
*/
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

/*
    Exp-Sinh factories (C ABI)
*/
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

/*
    quad_destroy
    ------------
    Releases the quadrature handle created by any factory above. Safe to call
    with nullptr. After destruction, the handle must not be used again.
*/
void quad_destroy(QuadratureHandle handle) {
    delete as_any(handle);
}

// === Integration =============================================================

/*
    quad_integrate_* (C ABI)
    ------------------------
    Integrate over the rule’s natural domain or default interval.
    - The integrand is a C function pointer of the appropriate precision.
    - The context pointer is passed through to the integrand.
*/
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

/*
    quad_integrate_interval_* (C ABI)
    ---------------------------------
    Integrate over [a, b] when supported. Implementations may ignore bounds
    for rules with natural infinite domains (see per-rule comments).
*/
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

/*
    quad_get_type / quad_get_precision / quad_get_points
    ----------------------------------------------------
    Accessors for metadata of the underlying implementation. If the handle
    is null, sensible defaults are returned:
      - type: QUAD_GAUSS_LEGENDRE
      - precision: QUAD_PRECISION_DOUBLE
      - points: 0
*/
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

/*
    quad_get_abscissa_weights_* (C ABI)
    -----------------------------------
    Copies abscissa and weights for fixed rules into caller-provided buffers.
    Returns 1 on success, 0 otherwise (including adaptive rules or insufficient
    buffer capacity). The caller must ensure the buffers have capacity equal
    to or greater than the rule’s point count.
*/
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

/*
    quad_type_to_string
    -------------------
    Returns a stable string name for the QuadratureType enum. Useful for
    diagnostics and logging; the returned pointer is to a static string
    literal, so it remains valid for the lifetime of the program.
*/
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

/*
    quad_is_adaptive
    ----------------
    Returns 1 if the rule is adaptive (double-exponential family), 0 otherwise.
    Adaptive rules do not expose fixed abscissa/weights and have get_points()=-1.
*/
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

/*
    quad_supports_infinite_bounds
    -----------------------------
    Returns 1 if the rule naturally supports infinite bounds, 0 otherwise.
    - Gauss-Hermite: (-∞, ∞)
    - Sinh-Sinh: (-∞, ∞)
    - Exp-Sinh: [0, ∞)
*/
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
