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
    /// An Exponential probability distribution.
    ///
    /// Definition and properties (rate form):
    /// - Parameter:
    ///   - λ (lambda): The rate parameter, strictly positive.
    /// - Support: [0, +∞).
    /// - PDF: f(x) = λ e^(−λ x) for x ≥ 0; 0 otherwise.
    /// - CDF: F(x) = 1 − e^(−λ x) for x ≥ 0.
    /// - Survival (SF): S(x) = e^(−λ x) for x ≥ 0.
    /// - Quantile: Q(p) = −ln(1 − p) / λ for p in [0, 1).
    /// - Hazard: h(x) = λ (constant, memoryless property).
    /// - Cumulative hazard: H(x) = λ x.
    ///
    /// Summary statistics (for λ > 0):
    /// - Mean: 1 / λ
    /// - Variance: 1 / λ²
    /// - Mode: 0
    /// - Median: ln 2 / λ
    /// - Skewness: 2
    /// - Kurtosis (Pearson): 9 (excess = 6)
    /// - Differential entropy (nats): 1 − ln λ
    ///
    /// Scale form:
    /// - Alternatively parameterized by a scale θ = 1 / λ (θ > 0), with PDF f(x) = (1/θ) e^(−x/θ).
    ///
    /// Implementation notes:
    /// - This is a thin, value-type Swift wrapper delegating core computations to a dynamic backend
    ///   represented by ``Distribution.Dynamic``. All stored properties are immutable; the type is
    ///   `Sendable` and safe to pass across concurrency domains.
    ///
    /// - Generic parameter `T`: The floating-point type (`Float`, `Double`, or platform-available `Float80`).
    ///
    /// Errors:
    /// - Throws ``DistributionError/parameterNotPositive(name:value:)`` when a positive parameter
    ///   requirement is violated (e.g., `lambda ≤ 0` or `scale ≤ 0`) during initialization.
    ///
    /// See also:
    /// - Wikipedia: Exponential distribution
    /// - Boost.Math: exponential_distribution
    public struct Exponential<T: Real & BinaryFloatingPoint & Sendable>: Sendable, DistributionProtocol {
        /// The real number type used by this distribution.
        typealias RealType = T

        /// The rate parameter λ, strictly positive.
        ///
        /// - Interpretation: Expected value is `1 / λ` and the hazard is the constant `λ`.
        public let lambda: T

        /// Internal dynamic backend that provides numerically robust implementations
        /// of the distribution’s core operations (PDF/CDF/SF, quantiles, hazards, and moments).
        private let dyn: Distribution.Dynamic<T>

        /// Creates an Exponential distribution with the given rate parameter λ.
        ///
        /// - Parameter lambda: The positive rate parameter (λ > 0).
        ///
        /// - Throws:
        ///   - ``DistributionError/parameterNotPositive(name:value:)`` if `lambda ≤ 0`.
        ///
        /// - Examples:
        ///   - `let exp1 = try Distribution.Exponential<Double>(lambda: 2.0)`
        ///   - `let x = try exp1.quantile(0.95)`
        public init(lambda: T) throws {
            guard lambda > 0 else {
                throw DistributionError.parameterNotPositive(
                    name: "lambda",
                    value: lambda
                )
            }
            self.lambda = lambda
            self.dyn = try Distribution.Dynamic<T>(
                distributionName: "exponential",
                parameters: [
                    "lambda": lambda
                ]
            )
        }

        /// Creates an Exponential distribution with the given scale parameter θ.
        ///
        /// - Parameter scale: The positive scale parameter θ = 1 / λ (θ > 0).
        ///
        /// - Throws:
        ///   - ``DistributionError/parameterNotPositive(name:value:)`` if `scale ≤ 0`.
        ///
        /// - Example:
        ///   ```swift
        ///   // θ = 0.5 ⇒ λ = 2
        ///   let expScale = try Distribution.Exponential<Double>(scale: 0.5)
        ///   ```
        public init(scale: T) throws {
            guard scale > 0 else {
                throw DistributionError.parameterNotPositive(
                    name: "scale",
                    value: scale
                )
            }
            self.lambda = 1 / scale
            self.dyn = try Distribution.Dynamic<T>(
                distributionName: "exponential",
                parameters: [
                    "lambda": 1 / scale
                ]
            )
        }
        
        // MARK: DistributionProtocol — Support

        /// The lower bound of the support.
        ///
        /// For Exponential, this is `0`.
        public var supportLowerBound: T { dyn.supportLowerBound }

        /// The upper bound of the support.
        ///
        /// For Exponential, this is `+∞`.
        public var supportUpperBound: T { dyn.supportUpperBound }

        /// The support range as a tuple `(lower, upper)`.
        ///
        /// For Exponential, this is typically `[0, +∞)`.
        public var range: (lower: T, upper: T) { dyn.range }

        // MARK: PDF/CDF/SF/Quantile

        /// The probability density function (PDF) evaluated at `x`.
        ///
        /// - Parameter x: Evaluation point (usually `x ≥ 0`; the density is `0` for `x < 0`).
        /// - Returns: The density value `f(x)`.
        /// - Throws: Rethrows any underlying backend error.
        public func pdf(_ x: T) throws -> T { try dyn.pdf(x) }

        /// The natural logarithm of the PDF evaluated at `x`.
        ///
        /// - Parameter x: Evaluation point.
        /// - Returns: The log-density `log f(x)`.
        /// - Throws: Rethrows any underlying backend error.
        public func logPdf(_ x: T) throws -> T { try dyn.logPdf(x) }

        /// The cumulative distribution function (CDF) evaluated at `x`.
        ///
        /// - Parameter x: Evaluation point.
        /// - Returns: `P(X ≤ x)`.
        /// - Throws: Rethrows any underlying backend error.
        public func cdf(_ x: T) throws -> T { try dyn.cdf(x) }

        /// The survival function (SF) evaluated at `x`.
        ///
        /// - Parameter x: Evaluation point.
        /// - Returns: `P(X > x) = 1 − CDF(x)`.
        /// - Throws: Rethrows any underlying backend error.
        public func sf(_ x: T) throws -> T { try dyn.sf(x) }

        /// The lower-tail quantile function (inverse CDF) evaluated at probability `p`.
        ///
        /// - Parameter p: A probability in `[0, 1]`.
        /// - Returns: The value `x` such that `P(X ≤ x) = p`.
        /// - Throws: If `p` is outside `[0, 1]` or if the computation fails.
        public func quantile(_ p: T) throws -> T { try dyn.quantile(p) }

        /// The upper-tail quantile function (inverse SF) evaluated at probability `q`.
        ///
        /// - Parameter q: A probability in `[0, 1]` representing the upper tail probability.
        /// - Returns: The value `x` such that `P(X > x) = q`.
        /// - Throws: If `q` is outside `[0, 1]` or if the computation fails.
        public func quantileComplement(_ q: T) throws -> T {
            try dyn.quantileComplement(q)
        }

        // MARK: Moments

        /// The mean of the distribution, if defined.
        ///
        /// - Note: For Exponential, `E[X] = 1 / λ` (always defined for `λ > 0`).
        public var mean: T? { dyn.mean }

        /// The variance of the distribution, if defined.
        ///
        /// - Note: For Exponential, `Var[X] = 1 / λ²`.
        public var variance: T? { dyn.variance }

        /// The mode of the distribution, if defined.
        ///
        /// - Note: For Exponential, `mode = 0`.
        public var mode: T? { dyn.mode }

        /// The median of the distribution.
        ///
        /// - Note: For Exponential, `median = ln 2 / λ`.
        public var median: T { dyn.median }

        /// The skewness of the distribution, if defined.
        ///
        /// - Note: For Exponential, `skewness = 2`.
        public var skewness: T? { dyn.skewness }

        /// The kurtosis of the distribution (Pearson), if defined.
        ///
        /// - Note: For Exponential, `kurtosis = 9` (excess = 6).
        public var kurtosis: T? { dyn.kurtosis }

        /// The excess kurtosis of the distribution, if defined.
        ///
        /// - Note: For Exponential, `excess = 6`.
        public var kurtosisExcess: T? { dyn.kurtosisExcess }

        // MARK: Hazards

        /// The (instantaneous) hazard function `h(x) = f(x) / S(x)` evaluated at `x`.
        ///
        /// - Parameter x: Evaluation point.
        /// - Returns: The hazard at `x`. For Exponential, this is the constant `λ`.
        /// - Throws: Rethrows any underlying backend error.
        public func hazard(_ x: T) throws -> T { try dyn.hazard(x) }

        /// The cumulative hazard function `H(x) = −log S(x)` evaluated at `x`.
        ///
        /// - Parameter x: Evaluation point.
        /// - Returns: The cumulative hazard at `x`. For Exponential, this equals `λ x`.
        /// - Throws: Rethrows any underlying backend error.
        public func chf(_ x: T) throws -> T { try dyn.chf(x) }

        // Lattice/discrete-only properties (continuous ⇒ nil)

        /// For lattice (discrete) distributions, the spacing between adjacent support points.
        ///
        /// - Note: Exponential is continuous; this is `nil`.
        public var latticeStep: T? { nil }

        /// For lattice (discrete) distributions, the origin (offset) of the lattice.
        ///
        /// - Note: Exponential is continuous; this is `nil`.
        public var latticeOrigin: T? { nil }

        /// The (differential) entropy h[X] where defined.
        ///
        /// - Note: For Exponential with rate λ, `h = 1 − ln λ` (nats). The backend provides this value.
        public var entropy: T? { dyn.entropy }
    }
}
