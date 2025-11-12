//
//  Created by Volker Thieme
//  Copyright © 2025 Volker Thieme. All rights reserved.
//
//  See the project root LICENSE for licensing information.
//

import SwiftyBoostPrelude

extension Distribution {
    /// Scaled inverse chi-squared probability distribution.
    ///
    /// Let `ν > 0` denote the degrees of freedom and `τ > 0` the scale parameter.
    /// The density is supported on `x ∈ (0, +∞)` and is given by
    ///
    /// ```
    /// f(x; ν, τ) = ((ν τ / 2)^(ν/2) / Γ(ν / 2)) · x^(−ν/2 − 1) · exp(−(ν τ)/(2x)).
    /// ```
    ///
    /// The type is a thin Swift value wrapper over Boost.Math's
    /// ``boost::math::inverse_chi_squared_distribution`` with complete coverage
    /// of the standard distribution protocol (PDF/CDF/quantiles, hazards, summary
    /// statistics, and KL divergence). All computations are delegated to the dynamic
    /// distribution factory, ensuring consistent numerical behaviour across precisions.
    public struct InverseChiSquared<T: Real & BinaryFloatingPoint & Sendable>: Sendable, DistributionProtocol {
        /// Concrete floating-point type used by the distribution.
        public typealias RealType = T

        /// Degrees of freedom `ν`.
        public let degreesOfFreedom: T

        /// Positive scale parameter `τ`.
        public let scale: T

        private let dyn: Distribution.Dynamic<T>

        /// Creates an inverse chi-squared distribution.
        ///
        /// - Parameters:
        ///   - degreesOfFreedom: The positive degrees of freedom `ν`.
        ///   - scale: Optional positive scale parameter `τ`. When `nil`, the canonical
        ///     choice `1 / ν` is used.
        ///
        /// - Throws:
        ///   - ``DistributionError/parameterNotPositive(name:value:)`` if `degreesOfFreedom ≤ 0`
        ///     or the provided `scale ≤ 0`.
        ///   - ``DistributionError/parameterNotFinite(name:value:)`` if any argument is `NaN`
        ///     or infinite.
        public init(degreesOfFreedom: T, scale: T? = nil) throws {
            guard degreesOfFreedom.isFinite else {
                throw DistributionError<T>.parameterNotFinite(
                    name: "degreesOfFreedom",
                    value: degreesOfFreedom
                )
            }
            guard degreesOfFreedom > 0 else {
                throw DistributionError<T>.parameterNotPositive(
                    name: "degreesOfFreedom",
                    value: degreesOfFreedom
                )
            }

            let resolvedScale: T
            if let explicitScale = scale {
                guard explicitScale.isFinite else {
                    throw DistributionError<T>.parameterNotFinite(
                        name: "scale",
                        value: explicitScale
                    )
                }
                guard explicitScale > 0 else {
                    throw DistributionError<T>.parameterNotPositive(
                        name: "scale",
                        value: explicitScale
                    )
                }
                resolvedScale = explicitScale
            } else {
                resolvedScale = 1 / degreesOfFreedom
            }

            self.degreesOfFreedom = degreesOfFreedom
            self.scale = resolvedScale

            let parameters: [String: T] = [
                "df": degreesOfFreedom,
                "scale": resolvedScale
            ]
            self.dyn = try Distribution.Dynamic<T>(
                distributionName: "inverse_chi_squared",
                parameters: parameters
            )
        }

        // MARK: - Support

        /// Lower bound of the support (`0` for the inverse chi-squared distribution).
        public var supportLowerBound: T { dyn.supportLowerBound }

        /// Upper bound of the support (`+∞` in exact arithmetic).
        public var supportUpperBound: T { dyn.supportUpperBound }

        /// Convenience tuple containing both support bounds.
        public var range: (lower: T, upper: T) { dyn.range }

        // MARK: - Core functions

        /// Probability density function `f(x)` evaluated at `x`.
        public func pdf(_ x: T) throws -> T { try dyn.pdf(x) }

        /// Natural logarithm of the PDF evaluated at `x`.
        public func logPdf(_ x: T) throws -> T { try dyn.logPdf(x) }

        /// Cumulative distribution function `F(x) = P(X ≤ x)`.
        public func cdf(_ x: T) throws -> T { try dyn.cdf(x) }

        /// Survival function `S(x) = P(X > x)`.
        public func sf(_ x: T) throws -> T { try dyn.sf(x) }

        /// Lower-tail quantile (inverse CDF) at probability `p`.
        public func quantile(_ p: T) throws -> T { try dyn.quantile(p) }

        /// Upper-tail quantile (inverse survival function) at probability `q`.
        public func quantileComplement(_ q: T) throws -> T {
            try dyn.quantileComplement(q)
        }

        // MARK: - Moments

        /// Mean `E[X]`, where defined (`nil` for `ν ≤ 2`).
        public var mean: T? { dyn.mean }

        /// Variance `Var[X]`, where defined (`nil` for `ν ≤ 4`).
        public var variance: T? { dyn.variance }

        /// Mode of the distribution, when defined.
        public var mode: T? { dyn.mode }

        /// Median of the distribution.
        public var median: T { dyn.median }

        /// Skewness, when defined (`nil` for `ν ≤ 6`).
        public var skewness: T? { dyn.skewness }

        /// Kurtosis (Pearson’s β₂), when defined (`nil` for `ν ≤ 8`).
        public var kurtosis: T? { dyn.kurtosis }

        /// Excess kurtosis (β₂ − 3), when defined.
        public var kurtosisExcess: T? { dyn.kurtosisExcess }

        // MARK: - Hazards

        /// Instantaneous hazard function `h(x) = f(x) / S(x)`.
        public func hazard(_ x: T) throws -> T { try dyn.hazard(x) }

        /// Cumulative hazard function `H(x) = −log S(x)`.
        public func chf(_ x: T) throws -> T { try dyn.chf(x) }

        // MARK: - Lattice metadata

        /// Lattice spacing for discrete distributions (not applicable here).
        public var latticeStep: T? { nil }

        /// Lattice origin for discrete distributions (not applicable here).
        public var latticeOrigin: T? { nil }

        // MARK: - Entropy & divergence

        /// Differential entropy `h[X]` (nats), when available.
        ///
        /// Boost.Math does not currently expose this quantity; the value is `nil`.
        public var entropy: T? { dyn.entropy }

        /// Indicates whether the distribution is discrete (`true`) or continuous (`false`).
        ///
        /// The inverse chi-squared distribution is continuous, so this is always `false`.
        public var isDiscrete: Bool { dyn.isDiscrete }

        /// Computes the Kullback–Leibler divergence `D_KL(self || other)`.
        ///
        /// - Parameters:
        ///   - other: Reference inverse chi-squared distribution *Q*.
        ///   - options: Numerical integration configuration. Defaults to
        ///     ``Distribution/KLDivergenceOptions/automatic()`` for the concrete type `T`.
        /// - Returns: The divergence in nats, `nil` when undefined, or `Double.infinity`
        ///   when the integral diverges.
        public func klDivergence<D>(
            relativeTo other: D,
            options: Distribution.KLDivergenceOptions<T>
        ) throws -> T? where D: DistributionProtocol, D.RealType == T {
            try DistributionKLDivergenceHelper.evaluate(lhs: self, rhs: other, options: options)
        }
}
}
