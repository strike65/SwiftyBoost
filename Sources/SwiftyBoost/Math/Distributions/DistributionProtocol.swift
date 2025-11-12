//
//  Created by Volker Thieme
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

/// A common interface for real-valued probability distributions.
///
/// Conforming types provide:
/// - The distribution's support (domain where the density/mass is nonzero).
/// - Core functions such as the PDF, log-PDF, CDF, survival function, hazard, and cumulative hazard.
/// - Inverse functions (quantiles) for both lower- and upper-tail probabilities.
/// - Summary statistics such as mean, variance, median, skewness, and kurtosis.
///
/// Unless otherwise noted, methods are expected to operate on values within the
/// distribution's support. Calling them with values outside the support may throw
/// an error specific to the conforming type (for example, a domain error).

import SwiftyBoostPrelude

extension Distribution {
    /// Configuration for numerically estimating KL divergence when no closed form is available.
    public struct KLDivergenceOptions<T: Real & BinaryFloatingPoint & Sendable>: Sendable, Equatable {
        public var finiteRule: Quadrature.Rule
        public var semiInfiniteRule: Quadrature.Rule
        public var infiniteRule: Quadrature.Rule
        public var densityFloor: T
        public var discreteTailCutoff: T
        public var maxDiscreteEvaluations: Int
        /// Optional override for the lower integration bound. When set, continuous integrals and discrete sums
        /// start at `max(sharedSupportLower, integrationLowerBound)`.
        public var integrationLowerBound: T?
        /// Optional override for the upper integration bound. When set, continuous integrals and discrete sums
        /// stop at `min(sharedSupportUpper, integrationUpperBound)`.
        public var integrationUpperBound: T?

        /// Creates a configuration for KL divergence evaluation.
        ///
        /// - Parameters:
        ///   - finiteRule: Quadrature rule for finite intervals.
        ///   - semiInfiniteRule: Quadrature rule for semi-infinite tails.
        ///   - infiniteRule: Quadrature rule for integrals spanning the full real line.
        ///   - densityFloor: Minimum density used to avoid division by zero.
        ///   - discreteTailCutoff: Tail probability threshold for discrete truncation.
        ///   - maxDiscreteEvaluations: Safety cap for discrete summations.
        ///   - integrationLowerBound: Optional override for the lower integration limit.
        ///   - integrationUpperBound: Optional override for the upper integration limit.
        public init(
            finiteRule: Quadrature.Rule = .gaussKronrod(points: 61),
            semiInfiniteRule: Quadrature.Rule = .expSinh(),
            infiniteRule: Quadrature.Rule = .tanhSinh(),
            densityFloor: T? = nil,
            discreteTailCutoff: T? = nil,
            maxDiscreteEvaluations: Int = 250_000,
            integrationLowerBound: T? = nil,
            integrationUpperBound: T? = nil
        ) {
            self.finiteRule = finiteRule
            self.semiInfiniteRule = semiInfiniteRule
            self.infiniteRule = infiniteRule
            self.densityFloor = densityFloor ?? Self.defaultDensityFloor(for: T.self)
            self.discreteTailCutoff = discreteTailCutoff ?? Self.defaultTailCutoff(for: T.self)
            self.maxDiscreteEvaluations = maxDiscreteEvaluations
            self.integrationLowerBound = integrationLowerBound
            self.integrationUpperBound = integrationUpperBound
        }

        public static func automatic(for type: T.Type = T.self) -> Self {
            Self(
                finiteRule: .gaussKronrod(points: 61),
                semiInfiniteRule: .expSinh(),
                infiniteRule: .tanhSinh(),
                densityFloor: defaultDensityFloor(for: type),
                discreteTailCutoff: defaultTailCutoff(for: type),
                integrationLowerBound: nil,
                integrationUpperBound: nil
            )
        }

        private static func defaultDensityFloor(for type: T.Type) -> T {
            if type == Float.self {
                return T(1e-9)
            }
            if type == Double.self {
                return T(1e-18)
            }
            return T(1e-24)
        }

        private static func defaultTailCutoff(for type: T.Type) -> T {
            if type == Float.self {
                return T(1e-6)
            }
            return T(1e-9)
        }
    }
}
public protocol DistributionProtocol: Sendable {
    /// The real number type used by this distribution (for example, `Float`, `Double`, `Float80`).
    ///
    /// Conforming types should choose a floating-point type that suits their precision
    /// and performance needs. The associated type must be `Sendable` so values can be
    /// safely passed across concurrency domains.
    associatedtype RealType: Real & BinaryFloatingPoint & Sendable

    // MARK: - Support

    /// The lower bound of the distribution's support (domain where the distribution is defined).
    ///
    /// Depending on the distribution, the lower bound may be finite or infinite (represented
    /// by `-infinity` if supported by [Real](https://github.com/apple/swift-numerics)). Whether the bound is inclusive or exclusive
    /// is distribution-specific and should be documented by the conforming type.
    var supportLowerBound: RealType { get }

    /// The upper bound of the distribution's support (domain where the distribution is defined).
    ///
    /// Depending on the distribution, the upper bound may be finite or infinite (represented
    /// by `+infinity` if supported by [Real](https://github.com/apple/swift-numerics)). Whether the bound is inclusive or exclusive
    /// is distribution-specific and should be documented by the conforming type.
    var supportUpperBound: RealType { get }

    // MARK: - Core functions

    /// The cumulative distribution function (CDF).
    ///
    /// For continuous distributions, this is F(x) = P(X ≤ x).
    ///
    /// - Parameter x: The evaluation point.
    /// - Returns: The CDF value in [0, 1].
    /// - Throws: If `x` is outside the supported domain or if the computation fails.
    func cdf(_ x: RealType) throws -> RealType

    /// The natural logarithm of the probability density function (log-PDF).
    ///
    /// Implementations should prefer numerically stable computations and avoid computing
    /// `log(pdf(x))` directly where possible.
    ///
    /// - Parameter x: The evaluation point.
    /// - Returns: The log-PDF at `x`.
    /// - Throws: If `x` is outside the supported domain or if the computation fails.
    func logPdf(_ x: RealType) throws -> RealType

    /// The probability density function (PDF).
    ///
    /// - Parameter x: The evaluation point.
    /// - Returns: The PDF at `x` (nonnegative, integrating to 1 over the support).
    /// - Throws: If `x` is outside the supported domain or if the computation fails.
    func pdf(_ x: RealType) throws -> RealType

    /// The survival function (SF), also known as the complementary CDF.
    ///
    /// For continuous distributions, this is S(x) = P(X > x). Where numerically advantageous,
    /// implementations should compute this directly rather than using `1 - cdf(x)`.
    ///
    /// - Parameter x: The evaluation point.
    /// - Returns: The survival probability in [0, 1].
    /// - Throws: If `x` is outside the supported domain or if the computation fails.
    func sf(_ x: RealType) throws -> RealType

    /// The (instantaneous) hazard function.
    ///
    /// For continuous distributions, h(x) = f(x) / S(x) where defined.
    ///
    /// - Parameter x: The evaluation point.
    /// - Returns: The hazard at `x`.
    /// - Throws: If `x` is outside the supported domain or if the computation fails.
    func hazard(_ x: RealType) throws -> RealType

    /// The cumulative hazard function (CHF).
    ///
    /// For continuous distributions, H(x) = −log S(x) where defined.
    ///
    /// - Parameter x: The evaluation point.
    /// - Returns: The cumulative hazard at `x`.
    /// - Throws: If `x` is outside the supported domain or if the computation fails.
    func chf(_ x: RealType) throws -> RealType

    // MARK: - Inverses

    /// The lower-tail quantile function, i.e. the inverse CDF.
    ///
    /// Returns the value `x` such that F(x) = p for `p` in [0, 1].
    /// Implementations should be monotonic in `p` and numerically stable near 0 and 1.
    ///
    /// - Parameter p: A probability in [0, 1].
    /// - Returns: The quantile corresponding to probability `p`.
    /// - Throws: If `p` is outside [0, 1] or if the computation fails to converge.
    func quantile(_ p: RealType) throws -> RealType

    /// The upper-tail quantile function, i.e. the inverse survival function.
    ///
    /// Returns the value `x` such that S(x) = p for `p` in [0, 1].
    /// Implementations should be monotonic in `p` and numerically stable near 0 and 1.
    ///
    /// - Parameter p: An upper-tail probability in [0, 1].
    /// - Returns: The upper-tail quantile corresponding to probability `p`.
    /// - Throws: If `p` is outside [0, 1] or if the computation fails to converge.
    func quantileComplement(_ p: RealType) throws -> RealType

    // MARK: - Summary statistics

    /// The expected value (mean) of the distribution.
    ///
    /// If the mean does not exist (diverges), conforming types should document and choose
    /// an appropriate behavior (for example, `±infinity`, `nan`, or `nil`).
    var mean: RealType? { get }

    /// The mode (value at which the PDF attains its maximum), if it exists and is unique.
    ///
    /// Some distributions do not have a unique mode; in such cases, this may be `nil`.
    var mode: RealType? { get }

    /// The variance of the distribution.
    ///
    /// If the variance does not exist (diverges), conforming types should document and choose
    /// an appropriate behavior (for example, `±infinity`, `nan`, or `nil`).
    var variance: RealType? { get }

    /// The median of the distribution.
    ///
    /// For continuous unimodal distributions this is typically unique; for others,
    /// a conventional choice should be documented by the conforming type.
    var median: RealType { get }

    /// The skewness of the distribution.
    ///
    /// Defined as E[((X − μ)/σ)^3] where it exists.
    var skewness: RealType? { get }

    /// The kurtosis of the distribution (Pearson's definition).
    ///
    /// Defined as E[((X − μ)/σ)^4] where it exists. This value equals 3 for the normal distribution.
    var kurtosis: RealType? { get }

    /// The excess kurtosis of the distribution.
    ///
    /// Defined as `kurtosis − 3`. Equals 0 for the normal distribution.
    var kurtosisExcess: RealType? { get }

    /// A convenient representation of the distribution's overall range.
    ///
    /// Unless documented otherwise by the conforming type, this is typically equivalent
    /// to `(supportLowerBound, supportUpperBound)`, which may include infinite values.
    var range: (lower: RealType, upper: RealType) { get }

    /// For lattice (discrete) distributions, the spacing between adjacent support points.
    ///
    /// Continuous distributions should return `nil`.
    var latticeStep: RealType? { get }

    /// For lattice (discrete) distributions, the origin (offset) of the lattice.
    ///
    /// Continuous distributions should return `nil`.
    var latticeOrigin: RealType? { get }

    /// The entropy H[X] (Shannon) where defined.
    ///
    /// Implementations should return `nil` when the entropy is undefined or does not
    /// exist (for example, diverges). Units are in nats (natural logarithm base).
    var entropy: RealType? { get }

    /// The Kullback–Leibler divergence `D_KL(self || other)` where defined.
    ///
    /// - Parameters:
    ///   - other: Distribution to compare against (must share support).
    ///   - options: Numerical integration/summation configuration.
    /// - Returns: Divergence in nats, `nil` when undefined, or `infinity` when divergent.
    func klDivergence<D: DistributionProtocol>(
        relativeTo other: D,
        options: Distribution.KLDivergenceOptions<RealType>
    ) throws -> RealType? where D.RealType == RealType

    /// Indicates whether the distribution is discrete (`true`) or continuous (`false`).
    ///
    /// Implementations should return `true` when the support is lattice-based and `false` otherwise.
    var isDiscrete: Bool { get }
}

public extension DistributionProtocol {
    /// Default configuration derived from the distribution’s support.
    ///
    /// - Returns: ``Distribution/KLDivergenceOptions`` with integration limits aligned to the
    ///   finite support endpoints (when available) and the standard quadrature rules for the real line.
    func defaultKLDivergenceOptions() -> Distribution.KLDivergenceOptions<RealType> {
        var options = Distribution.KLDivergenceOptions<RealType>.automatic(for: RealType.self)
        let lower = supportLowerBound
        if lower.isFinite {
            options.integrationLowerBound = lower
        }
        let upper = supportUpperBound
        if upper.isFinite {
            options.integrationUpperBound = upper
        }
        return options
    }

    func klDivergence<D: DistributionProtocol>(
        relativeTo other: D,
        options: Distribution.KLDivergenceOptions<RealType>
    ) throws -> RealType? where D.RealType == RealType {
        try DistributionKLDivergenceHelper.evaluate(lhs: self, rhs: other, options: options)
    }

    /// Convenience overload that uses ``Distribution/KLDivergenceOptions/automatic(for:)``.
    func klDivergence<D: DistributionProtocol>(relativeTo other: D) throws -> RealType?
    where D.RealType == RealType {
        try klDivergence(relativeTo: other, options: defaultKLDivergenceOptions())
    }
}

@usableFromInline
enum DistributionKLDivergenceHelper {
    @usableFromInline
    static func evaluate<LHS: DistributionProtocol, RHS: DistributionProtocol>(
        lhs: LHS,
        rhs: RHS,
        options: Distribution.KLDivergenceOptions<LHS.RealType>
    ) throws -> LHS.RealType? where LHS.RealType == RHS.RealType {
        guard lhs.isDiscrete == rhs.isDiscrete else { return nil }
        if lhs.isDiscrete {
            return try discrete(lhs: lhs, rhs: rhs, options: options)
        } else {
            return try continuous(lhs: lhs, rhs: rhs, options: options)
        }
    }

    @usableFromInline
    static func continuous<LHS: DistributionProtocol, RHS: DistributionProtocol>(
        lhs: LHS,
        rhs: RHS,
        options: Distribution.KLDivergenceOptions<LHS.RealType>
    ) throws -> LHS.RealType? where LHS.RealType == RHS.RealType {
        typealias T = LHS.RealType
        guard let bounds = sharedBounds(
            lhs: lhs,
            rhs: rhs,
            lowerOverride: options.integrationLowerBound,
            upperOverride: options.integrationUpperBound
        ) else { return nil }
        if !(bounds.lower < bounds.upper) { return nil }
        let densityFloor = max(options.densityFloor, T.leastNonzeroMagnitude)
        let flag = DivergenceFlag<T>()
        let integrand: @Sendable (T) -> T = { x in
            let p = positiveDensity(lhs, at: x)
            if !p.isFinite || p <= densityFloor { return .zero }
            var q = positiveDensity(rhs, at: x)
            if !q.isFinite {
                flag.isInfinite = true
                return .zero
            }
            if q <= .zero {
                q = densityFloor
            }
            if q < densityFloor {
                q = densityFloor
            }
            let ratio = p / q
            if !ratio.isFinite || ratio <= .zero {
                return .zero
            }
            return p * T.log(ratio)
        }

        let value: T
        switch (bounds.lower.isFinite, bounds.upper.isFinite) {
        case (true, true):
            let lowerD = Double(bounds.lower)
            let upperD = Double(bounds.upper)
            guard lowerD.isFinite, upperD.isFinite, lowerD < upperD else { return nil }
            let integrator = try Quadrature.Integrator<T>(rule: options.finiteRule)
            value = try integrator
                .integrate(over: .finite(lower: lowerD, upper: upperD), integrand: integrand)
                .value
        case (true, false):
            let lower = bounds.lower
            guard lower.isFinite else { return nil }
            let integrator = try Quadrature.Integrator<T>(rule: options.semiInfiniteRule)
            value = try integrator.integrate { shift in
                integrand(lower + shift)
            }.value
        case (false, true):
            let upper = bounds.upper
            guard upper.isFinite else { return nil }
            let integrator = try Quadrature.Integrator<T>(rule: options.semiInfiniteRule)
            value = try integrator.integrate { shift in
                integrand(upper - shift)
            }.value
        default:
            let integrator = try Quadrature.Integrator<T>(rule: options.infiniteRule)
            value = try integrator.integrate(over: .automatic, integrand: integrand).value
        }

        if flag.isInfinite {
            return T.infinity
        }
        return value.isFinite ? value : nil
    }

    @usableFromInline
    static func discrete<LHS: DistributionProtocol, RHS: DistributionProtocol>(
        lhs: LHS,
        rhs: RHS,
        options: Distribution.KLDivergenceOptions<LHS.RealType>
    ) throws -> LHS.RealType? where LHS.RealType == RHS.RealType {
        typealias T = LHS.RealType
        guard let bounds = sharedBounds(
            lhs: lhs,
            rhs: rhs,
            lowerOverride: options.integrationLowerBound,
            upperOverride: options.integrationUpperBound
        ) else { return nil }
        let densityFloor = max(options.densityFloor, T.leastNonzeroMagnitude)
        let tailTolerance = max(options.discreteTailCutoff, T.leastNonzeroMagnitude)

        let startIndex: Int
        if bounds.lower.isInfinite {
            return nil
        } else if let s = discreteCeil(bounds.lower) {
            startIndex = s
        } else {
            return nil
        }

        var endIndex: Int? = nil
        if bounds.upper.isFinite {
            endIndex = discreteFloor(bounds.upper)
            if let end = endIndex, end < startIndex {
                return nil
            }
        }

        var idx = startIndex
        var divergence = T.zero
        var iterations = 0
        var infiniteFlag = false
        while true {
            if let end = endIndex, idx > end { break }
            let point = T(idx)
            let p = positiveDensity(lhs, at: point)
            if p > densityFloor {
                var q = positiveDensity(rhs, at: point)
                if !q.isFinite || q <= .zero {
                    infiniteFlag = true
                } else {
                    if q < densityFloor { q = densityFloor }
                    let ratio = p / q
                    if ratio.isFinite && ratio > .zero {
                        divergence += p * T.log(ratio)
                    }
                }
            }

            iterations += 1
            if iterations >= options.maxDiscreteEvaluations {
                return nil
            }

            if endIndex == nil {
                let tailSelf = positiveSurvival(lhs, above: point)
                let tailOther = positiveSurvival(rhs, above: point)
                if tailSelf <= tailTolerance && tailOther <= tailTolerance {
                    break
                }
            }

            idx += 1
        }

        if infiniteFlag {
            return T.infinity
        }
        return divergence
    }

    @usableFromInline
    static func sharedBounds<LHS: DistributionProtocol, RHS: DistributionProtocol>(
        lhs: LHS,
        rhs: RHS,
        lowerOverride: LHS.RealType?,
        upperOverride: LHS.RealType?
    ) -> (lower: LHS.RealType, upper: LHS.RealType)? where LHS.RealType == RHS.RealType {
        var lower = Swift.max(lhs.supportLowerBound, rhs.supportLowerBound)
        var upper = Swift.min(lhs.supportUpperBound, rhs.supportUpperBound)
        if let overrideLower = lowerOverride {
            lower = Swift.max(lower, overrideLower)
        }
        if let overrideUpper = upperOverride {
            upper = Swift.min(upper, overrideUpper)
        }
        if lower.isNaN || upper.isNaN || lower > upper {
            return nil
        }
        return (lower, upper)
    }

    @usableFromInline
    static func positiveDensity<D: DistributionProtocol>(
        _ distribution: D,
        at x: D.RealType
    ) -> D.RealType {
        guard let value = try? distribution.pdf(x), value.isFinite, value > 0 else {
            return .zero
        }
        return value
    }

    @usableFromInline
    static func positiveSurvival<D: DistributionProtocol>(
        _ distribution: D,
        above x: D.RealType
    ) -> D.RealType {
        guard let value = try? distribution.sf(x), value.isFinite else {
            return D.RealType.infinity
        }
        return value < 0 ? D.RealType.infinity : value
    }

    @usableFromInline
    static func discreteCeil<T: Real & BinaryFloatingPoint>(_ value: T) -> Int? {
        let dv = Double(value)
        guard dv.isFinite else { return nil }
        if dv >= Double(Int.max) { return nil }
        if dv <= Double(Int.min) { return nil }
        let rounded = dv.rounded(.up)
        return Int(rounded)
    }

    @usableFromInline
    static func discreteFloor<T: Real & BinaryFloatingPoint>(_ value: T) -> Int? {
        let dv = Double(value)
        guard dv.isFinite else { return nil }
        if dv >= Double(Int.max) { return Int(Int.max) }
        if dv <= Double(Int.min) { return nil }
        let rounded = dv.rounded(.down)
        return Int(rounded)
    }

    @usableFromInline
    final class DivergenceFlag<T: Real & BinaryFloatingPoint>: @unchecked Sendable {
        var isInfinite: Bool = false
    }
}
