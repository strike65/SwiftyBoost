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
    /// Gamma distribution Γ(k, θ) with shape k > 0 and scale θ > 0.
    ///
    /// Definitions
    /// - PDF: f(x; k, θ) = x^(k−1) e^(−x/θ) / (Γ(k) θ^k), for x ≥ 0
    /// - CDF (lower tail): F(x) = P(k, x/θ)
    /// - SF (upper tail): S(x) = Q(k, x/θ)
    /// - Quantile: F(x) = p ⇒ x = θ · P⁻¹(k, p)
    ///
    /// Moments (where defined)
    /// - Mean = kθ
    /// - Variance = kθ²
    /// - Mode = (k − 1)θ for k ≥ 1; undefined for k < 1
    /// - Skewness = 2 / √k
    /// - Kurtosis (Pearson) = 3 + 6/k
    /// - Excess kurtosis = 6/k
    ///
    /// Notes
    /// - This type constructs the underlying Boost `gamma_distribution` once and
    ///   keeps a small opaque handle; all evaluations reuse this instance.
    /// - Throws on invalid parameters (k ≤ 0 or θ ≤ 0) or invalid arguments to
    ///   evaluators (e.g., x < 0 for PDF/CDF/SF).
    public struct Gamma<T: Real & BinaryFloatingPoint & Sendable>: Sendable, DistributionProtocol {
        public typealias RealType = T

        /// Shape parameter k (> 0).
        public let shape: T
        /// Scale parameter θ (> 0).
        public let scale: T

        private let dyn: Distribution.Dynamic<T>

        /// Initialize Γ(k, θ).
        ///
        /// - Parameters:
        ///   - shape: k > 0
        ///   - scale: θ > 0 (default 1)
        /// - Throws: `DistributionError.parameterNotPositive` if parameters are invalid.
        public init(shape: T, scale: T = 1) throws {
            guard shape > 0 else {
                throw DistributionError.parameterNotPositive(name: "shape", value: shape)
            }
            guard scale > 0 else {
                throw DistributionError.parameterNotPositive(name: "scale", value: scale)
            }
            self.shape = shape
            self.scale = scale
            self.dyn = try Distribution.Dynamic<T>(
                distributionName: "gamma",
                parameters: [
                    "shape": shape,
                    "scale": scale,
                ]
            )
        }

        // MARK: DistributionProtocol — Support

        /// The lower bound of the support. For Γ(k, θ), this is 0.
        public var supportLowerBound: T {
            dyn.supportLowerBound
        }

        /// The upper bound of the support. For Γ(k, θ), this is +∞.
        public var supportUpperBound: T {
            dyn.supportUpperBound
        }

        /// A convenience tuple of the distribution’s overall range.
        ///
        /// For Γ(k, θ) this is typically (0, +∞), possibly adjusted by Boost’s
        /// internal policies for numerical underflow/overflow reporting.
        public var range: (lower: T, upper: T) { dyn.range }

        // MARK: PDF/CDF/SF/Quantile

        /// Evaluate the probability density function f(x) at `x`.
        ///
        /// - Parameter x: Must satisfy x ≥ 0 for Γ(k, θ).
        /// - Returns: The density value at `x`.
        /// - Throws: `DistributionError.parameterOutOfRange` if `x < 0`.
        public func pdf(_ x: T) throws -> T {
            guard x >= 0 else {
                throw DistributionError.parameterOutOfRange(
                    name: "x",
                    min: 0.0,
                    max: Double.infinity
                )
            }
            return try dyn.pdf(x)
        }

        /// Evaluate the lower-tail CDF F(x) = P(X ≤ x) at `x`.
        ///
        /// - Parameter x: Must satisfy x ≥ 0 for Γ(k, θ).
        /// - Returns: The lower-tail probability in [0, 1].
        /// - Throws: `DistributionError.parameterOutOfRange` if `x < 0`.
        public func cdf(_ x: T) throws -> T {
            guard x >= 0 else {
                throw DistributionError.parameterOutOfRange(
                    name: "x",
                    min: 0.0,
                    max: Double.infinity
                )
            }
            if x == 0 { return 0 }
            return try dyn.cdf(x)
        }

        /// Evaluate the survival function S(x) = P(X > x) at `x`.
        ///
        /// - Parameter x: Must satisfy x ≥ 0 for Γ(k, θ).
        /// - Returns: The upper-tail probability in [0, 1].
        /// - Throws: `DistributionError.parameterOutOfRange` if `x < 0`.
        public func sf(_ x: T) throws -> T {
            guard x >= 0 else {
                throw DistributionError.parameterOutOfRange(
                    name: "x",
                    min: 0.0,
                    max: Double.infinity
                )
            }
            if x == 0 { return 1 }
            return try dyn.sf(x)
        }

        /// Quantile (inverse lower-tail CDF): returns x such that F(x) = p.
        ///
        /// - Parameter p: Probability in [0, 1].
        /// - Returns: The lower-tail quantile x.
        /// - Throws: `DistributionError.parameterOutOfRange` if `p ∉ [0, 1]`.
        public func quantile(_ p: T) throws -> T {
            guard p >= 0 && p <= 1 else {
                throw DistributionError.parameterOutOfRange(
                    name: "p",
                    min: 0.0,
                    max: 1.0
                )
            }
            if p == 0 { return 0 }
            if p == 1 { return T.infinity }
            return try dyn.quantile(p)
        }

        /// Upper-tail quantile (inverse SF): returns x such that S(x) = q.
        ///
        /// - Parameter q: Upper-tail probability in [0, 1].
        /// - Returns: The upper-tail quantile x.
        /// - Throws: `DistributionError.parameterOutOfRange` if `q ∉ [0, 1]`.
        public func quantileComplement(_ q: T) throws -> T {
            guard q >= 0 && q <= 1 else {
                throw DistributionError.parameterOutOfRange(
                    name: "q",
                    min: 0.0,
                    max: 1.0
                )
            }
            if q == 1 { return 0 }
            if q == 0 { return T.infinity }
            return try dyn.quantileComplement(q)
        }

        // MARK: Moments

        /// The mean E[X] = kθ.
        public var mean: T? { dyn.mean }

        /// The variance Var[X] = kθ².
        public var variance: T? { dyn.variance }

        /// The mode, if defined.
        ///
        /// For Γ(k, θ), the mode is (k − 1)θ when k ≥ 1; undefined otherwise.
        /// Returns `nil` when not finite or undefined.
        public var mode: T? { dyn.mode }

        /// The hazard function h(x) = f(x) / S(x) where defined.
        ///
        /// - Parameter x: Must satisfy x ≥ 0 for Γ(k, θ).
        /// - Returns: The instantaneous hazard at `x`.
        /// - Throws: `DistributionError.parameterOutOfRange` if `x < 0`.
        public func hazard(_ x: T) throws -> T {
            guard x >= 0 else {
                throw DistributionError.parameterOutOfRange(
                    name: "x",
                    min: 0.0,
                    max: Double.infinity
                )
            }
            return try dyn.hazard(x)
        }

        /// The cumulative hazard H(x) = −log S(x) where defined.
        ///
        /// - Parameter x: Must satisfy x ≥ 0 for Γ(k, θ).
        /// - Returns: The cumulative hazard at `x`.
        /// - Throws: `DistributionError.parameterOutOfRange` if `x < 0`.
        public func chf(_ x: T) throws -> T {
            guard x >= 0 else {
                throw DistributionError.parameterOutOfRange(
                    name: "x",
                    min: 0.0,
                    max: Double.infinity
                )
            }
            return try dyn.chf(x)
        }

        /// The natural logarithm of the density, log f(x).
        ///
        /// - Parameter x: Must satisfy x ≥ 0 for Γ(k, θ).
        /// - Returns: The log-PDF at `x`.
        /// - Throws: `DistributionError.parameterOutOfRange` if `x < 0`.
        public func logPdf(_ x: T) throws -> T {
            let p = try pdf(x)
            return T(log(Double(p)))
        }

        // MARK: DistributionProtocol — Additional moments

        /// The median of the distribution (numerically evaluated by Boost).
        public var median: T { dyn.median }

        /// The skewness γ₁ = E[((X − μ)/σ)^3], where defined.
        ///
        /// For Γ(k, θ), γ₁ = 2 / √k.
        public var skewness: T? { dyn.skewness }

        /// The kurtosis (Pearson’s β₂) = E[((X − μ)/σ)^4], where defined.
        ///
        /// For Γ(k, θ), β₂ = 3 + 6/k. Equals 3 for a normal distribution.
        public var kurtosis: T? { dyn.kurtosis }

        /// The excess kurtosis γ₂ = kurtosis − 3, where defined.
        ///
        /// For Γ(k, θ), γ₂ = 6/k. Equals 0 for a normal distribution.
        public var kurtosisExcess: T? { dyn.kurtosisExcess }

        /// Lattice step for discrete distributions (nil for Γ).
        public var latticeStep: T? { nil }

        /// Lattice origin for discrete distributions (nil for Γ).
        public var latticeOrigin: T? { nil }

        /// Entropy H[X] (Shannon), where defined.
        public var entropy: T? {
            do {
                let lg = try SpecialFunctions.logGamma(self.shape)
                let dg = try SpecialFunctions.digamma(self.shape)
                return self.shape + T.log(self.scale) + lg + (1 - self.shape) * dg
            }
            catch _ {
                return dyn.entropy
            }
        }

        /// Indicates whether this distribution is discrete (`true`) or continuous (`false`).
        ///
        /// Gamma distributions are continuous on `[0, ∞)`, so this always returns `false`.
        public var isDiscrete: Bool { dyn.isDiscrete }

        /// Computes the Kullback–Leibler divergence `D_KL(self || other)` when defined.
        ///
        /// - Parameters:
        ///   - other: The reference gamma distribution *Q*.
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
