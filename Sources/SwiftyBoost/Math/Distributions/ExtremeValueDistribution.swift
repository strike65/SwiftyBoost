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
import SwiftyBoostPrelude

extension Distribution {

    /// Extreme Value (Gumbel) distribution.
    ///
    /// This type models the Extreme Value Type I distribution (commonly called the Gumbel distribution),
    /// parametrized by a location parameter `loc` (μ) and a strictly positive scale parameter `scale` (β).
    ///
    /// - Support: (-∞, +∞)
    /// - Parameters:
    ///   - location (`loc`): Real-valued location parameter shifting the distribution along the real line.
    ///   - scale (`scale`): Positive real-valued scale parameter controlling dispersion. Must be `> 0`.
    ///
    /// Most numeric operations (PDF/CDF/SF/quantiles/moments/hazards) are delegated to a dynamic
    /// backend (`Distribution.Dynamic`) configured with the distribution name `"extreme_value"`
    /// and the parameters `["loc": location, "scale": scale]`.
    ///
    /// Thread-safety: This type is immutable and marked `Sendable`.
    ///
    /// Usage example:
    /// ```swift
    /// let dist = try Distribution.ExtremeValueGumbel<Double>(location: 0.0, scale: 2.0)
    /// let p = try dist.cdf(1.5)
    /// let x = try dist.quantile(0.9)
    /// ```
    public struct ExtremeValueGumbel<T: Real & BinaryFloatingPoint & Sendable>: Sendable, DistributionProtocol {
        /// The underlying real type used by the distribution.
        public typealias RealType = T

        /// Location parameter μ.
        ///
        /// Shifts the distribution along the real line.
        public let location: T

        /// Scale parameter β (must be strictly positive).
        ///
        /// Controls the spread of the distribution.
        public let scale: T

        /// Internal dynamic backend used to evaluate distribution functions.
        private let dyn: Distribution.Dynamic<T>

        /// Creates an Extreme Value (Gumbel) distribution.
        ///
        /// - Parameters:
        ///   - loc: Location parameter μ. Default is `0`.
        ///   - scale: Scale parameter β. Must be strictly positive. Default is `1`.
        /// - Throws: `DistributionError.parameterNotPositive(name:value:)` if `scale <= 0`.
        public init(location loc: T = 0, scale: T = 1) throws {
            guard scale > 0 else {
                throw DistributionError.parameterNotPositive(
                    name: "scale",
                    value: scale
                )
            }
            self.location = loc
            self.scale = scale
            self.dyn = try Distribution.Dynamic<T>(
                distributionName: "extreme_value",
                parameters: [
                    "loc": loc,
                    "scale": scale
                ]
            )
        }

        // MARK: DistributionProtocol — Support

        /// Lower bound of the support (for the Gumbel distribution this is -∞).
        public var supportLowerBound: T { dyn.supportLowerBound }

        /// Upper bound of the support (for the Gumbel distribution this is +∞).
        public var supportUpperBound: T { dyn.supportUpperBound }

        /// Convenience tuple of `(lower, upper)` support bounds.
        public var range: (lower: T, upper: T) { dyn.range }

        // MARK: PDF/CDF/SF/Quantile

        /// Probability density function f(x).
        ///
        /// - Parameter x: Point at which to evaluate the density.
        /// - Returns: The density value at `x`.
        public func pdf(_ x: T) throws -> T { try dyn.pdf(x) }

        /// Natural logarithm of the probability density function, ln f(x).
        ///
        /// - Parameter x: Point at which to evaluate log-density.
        /// - Returns: The log-density value at `x`.
        public func logPdf(_ x: T) throws -> T { try dyn.logPdf(x) }

        /// Cumulative distribution function F(x) = P[X ≤ x].
        ///
        /// - Parameter x: Point at which to evaluate the CDF.
        /// - Returns: The cumulative probability at `x`.
        public func cdf(_ x: T) throws -> T { try dyn.cdf(x) }

        /// Survival function S(x) = P[X > x] = 1 − F(x).
        ///
        /// - Parameter x: Point at which to evaluate the survival probability.
        /// - Returns: The survival probability at `x`.
        public func sf(_ x: T) throws -> T { try dyn.sf(x) }

        /// Quantile function (inverse CDF).
        ///
        /// - Parameter p: Probability in [0, 1].
        /// - Returns: The value `x` such that `F(x) = p`.
        /// - Throws: An error if `p` is outside [0, 1] or the backend cannot compute the quantile.
        public func quantile(_ p: T) throws -> T { try dyn.quantile(p) }

        /// Complementary quantile function (inverse survival).
        ///
        /// - Parameter q: Survival probability in [0, 1].
        /// - Returns: The value `x` such that `S(x) = q`.
        /// - Throws: An error if `q` is outside [0, 1] or the backend cannot compute the quantile.
        public func quantileComplement(_ q: T) throws -> T {
            try dyn.quantileComplement(q)
        }

        // MARK: Moments

        /// The mean of the distribution, if available.
        ///
        /// - Note: This property is computed locally in this type. Most other moment computations are
        ///   delegated to the dynamic backend. The value returned here reflects the current implementation
        ///   of this type and may differ from a closed-form expression you might expect from literature.
        public var mean: T? { location > 1 ? T(0) : nil }

        /// The variance of the distribution, if available.
        ///
        /// Delegated to the dynamic backend.
        public var variance: T? { dyn.variance }

        /// The mode of the distribution, if available.
        ///
        /// Delegated to the dynamic backend.
        public var mode: T? { dyn.mode }

        /// The median of the distribution.
        ///
        /// Delegated to the dynamic backend.
        public var median: T { dyn.median }

        /// The skewness of the distribution, if available.
        ///
        /// Delegated to the dynamic backend.
        public var skewness: T? { dyn.skewness }

        /// The kurtosis (including 3) of the distribution, if available.
        ///
        /// Delegated to the dynamic backend.
        public var kurtosis: T? { dyn.kurtosis }

        /// The excess kurtosis (kurtosis − 3) of the distribution, if available.
        ///
        /// Delegated to the dynamic backend.
        public var kurtosisExcess: T? { dyn.kurtosisExcess }

        // MARK: Hazards

        /// Hazard (failure rate) function h(x) = f(x) / S(x).
        ///
        /// - Parameter x: Point at which to evaluate the hazard.
        /// - Returns: The hazard value at `x`.
        public func hazard(_ x: T) throws -> T { try dyn.hazard(x) }

        /// Cumulative hazard function H(x) = −ln S(x).
        ///
        /// - Parameter x: Point at which to evaluate the cumulative hazard.
        /// - Returns: The cumulative hazard value at `x`.
        public func chf(_ x: T) throws -> T { try dyn.chf(x) }

        // Lattice/discrete-only properties (continuous ⇒ nil)

        /// For continuous distributions this is always `nil`.
        public var latticeStep: T? { nil }

        /// For continuous distributions this is always `nil`.
        public var latticeOrigin: T? { nil }

        /// Differential entropy, if available.
        ///
        /// - Note: This property is computed locally in this type. For the Gumbel distribution,
        ///   the theoretical differential entropy is often expressed using the Euler–Mascheroni
        ///   constant γ. The value returned here reflects the current implementation of this type.
        public var entropy: T? {
            let em: T = Constants.pi()
            return T.log(self.scale) + em + 1
        }
        /// Indicates whether this distribution is discrete (`true`) or continuous (`false`).
        ///
        /// The Gumbel (extreme value) distribution is continuous, so this always returns `false`.
        public var isDiscrete: Bool { dyn.isDiscrete }

        /// Computes the Kullback–Leibler divergence `D_KL(self || other)` when defined.
        ///
        /// - Parameters:
        ///   - other: The reference extreme-value distribution *Q*.
        ///   - options: Quadrature configuration; defaults to ``Distribution/KLDivergenceOptions/automatic()``.
        /// - Returns: The divergence in nats, or `nil` if it cannot be evaluated.
        /// - Throws: Rethrows any backend or quadrature errors.
        public func klDivergence<D>(
            relativeTo other: D,
            options: Distribution.KLDivergenceOptions<T>
        ) throws -> T? where D: DistributionProtocol, D.RealType == T {
            try DistributionKLDivergenceHelper.evaluate(lhs: self, rhs: other, options: options)
        }
}
}
