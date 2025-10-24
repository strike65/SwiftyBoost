//
//  bs_quadrature_removed.hxx
//  CBoostBridge
//
//  Archived implementations for Gauss–Hermite, Gauss–Laguerre, and
//  Gauss–Jacobi quadrature. These rely on Boost ≥ 1.78. The current vendored
//  Boost release does not provide these integrators, so the code is kept here
//  for future reference without being compiled.
//

/*
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
        auto lambda = [f, context](Real x) {
            return f ? f(x, context) : Real(0);
        };
        result.result = GaussType::integrate(lambda);
        return result;
    }

    Result integrate_interval(Integrand f, void* context, Real, Real) override {
        return integrate(f, context);
    }

    QuadratureType get_type() const override { return QUAD_GAUSS_HERMITE; }
    int get_points() const override { return static_cast<int>(Points); }
};
*/

/*
template <typename Real, unsigned Points>
class GaussLaguerreImpl final : public QuadratureBase<Real> {
    using Base = QuadratureBase<Real>;
    using Result = typename Base::Result;
    using Integrand = typename Base::Integrand;
    using GaussType = gauss_laguerre<Real, Points>;

public:
    explicit GaussLaguerreImpl(Real alpha = Real(0)) : alpha_(alpha) {}

    Result integrate(Integrand f, void* context) override {
        Result result{};
        result.converged = 1;
        result.iterations = 1;
        result.function_calls = static_cast<int>(Points);
        auto lambda = [f, context](Real x) {
            return f ? f(x, context) : Real(0);
        };
        if (alpha_ == Real(0)) {
            result.result = GaussType::integrate(lambda);
        } else {
            gauss_laguerre<Real, Points> quad(alpha_);
            result.result = quad.integrate(lambda);
        }
        return result;
    }

    Result integrate_interval(Integrand f, void* context, Real, Real) override {
        return integrate(f, context);
    }

    QuadratureType get_type() const override { return QUAD_GAUSS_LAGUERRE; }
    int get_points() const override { return static_cast<int>(Points); }

private:
    Real alpha_;
};
*/

/*
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
        auto lambda = [f, context](Real x) {
            return f ? f(x, context) : Real(0);
        };
        result.result = quad.integrate(lambda, a, b);
        return result;
    }

    QuadratureType get_type() const override { return QUAD_GAUSS_JACOBI; }
    int get_points() const override { return static_cast<int>(Points); }

private:
    Real alpha_;
    Real beta_;
};
*/

