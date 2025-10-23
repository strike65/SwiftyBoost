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
    /// A Cauchy (Lorentz) probability distribution.
    ///
    /// Definition and properties:
    /// - Parameters:
    ///   - location (x0): The central location (median and mode).
    ///   - scale (γ): The positive scale parameter (half-width at half-maximum).
    /// - Support: (-∞, +∞).
    /// - Density: f(x) = 1 / [πγ (1 + ((x − x0)/γ)^2)].
    /// - CDF: F(x) = 1/π arctan((x − x0)/γ) + 1/2.
    /// - Quantile: Q(p) = x0 + γ tan(π(p − 1/2)) for p in (0, 1).
    ///
    /// Statistical notes:
    /// - The mean and all moments of order ≥ 1 are undefined (do not exist).
    /// - The variance is undefined (does not exist).
    /// - The median equals the location parameter.
    /// - The mode equals the location parameter.
    /// - Differential entropy (where defined by convention) is h = ln(4πγ).
    ///
    /// Implementation notes:
    /// - This is a thin, value-type Swift wrapper that delegates core computations (PDF, CDF, SF, quantiles,
    ///   hazards, and derived moments) to a dynamic backend represented by ``Distribution.Dynamic``.
    /// - All stored properties are immutable; the type is `Sendable` and thus safe to pass across concurrency domains.
    ///
    /// - Generic parameter `T`: The floating-point type (`Float`, `Double`, or platform-available `Float80`).
    ///
    /// Errors:
    /// - Throws ``DistributionError/parameterNotPositive(name:value:)`` when `scale ≤ 0` during initialization.
    ///
    /// See also:
    /// - Wikipedia: Cauchy distribution
    /// - Boost.Math: cauchy_distribution
    public struct Cauchy<T: Real & BinaryFloatingPoint & Sendable>: Sendable, DistributionProtocol {
        /// The real number type used by this distribution.
        typealias RealType = T

        /// The location parameter x0 (any finite real).
        ///
        /// - Interpretation: Also equals the median and the mode of the distribution.
        public let location: T

        /// The positive scale parameter γ.
        ///
        /// - Requirement: `scale > 0`.
        /// - Interpretation: Half-width at half-maximum (HWHM).
        public let scale: T

        /// Internal dynamic backend that provides numerically robust implementations
        /// of the distribution’s core operations.
        private let dyn: Distribution.Dynamic<T>

        /// Creates a Cauchy distribution with the given location and scale.
        ///
        /// - Parameters:
        ///   - loc: The location parameter x0 (finite real). Default is `0`.
        ///   - scale: The positive scale parameter γ. Must satisfy `scale > 0`. Default is `1`.
        ///
        /// - Throws:
        ///   - ``DistributionError/parameterNotPositive(name:value:)`` if `scale ≤ 0`.
        ///
        /// - Examples:
        ///   - Standard Cauchy:
        ///     - `let dist = try Distribution.Cauchy<Double>()`
        ///   - Shifted/scaled:
        ///     - `let dist = try Distribution.Cauchy<Double>(location: 2.5, scale: 0.75)`
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
                distributionName: "cauchy",
                parameters: [
                    "loc": loc,
                    "scale": scale
                ]
            )
        }

        // MARK: DistributionProtocol — Support

        /// The lower bound of the support.
        ///
        /// For Cauchy, this is `-∞`.
        public var supportLowerBound: T { dyn.supportLowerBound }

        /// The upper bound of the support.
        ///
        /// For Cauchy, this is `+∞`.
        public var supportUpperBound: T { dyn.supportUpperBound }

        /// The support range as a tuple `(lower, upper)`.
        ///
        /// For Cauchy, this is typically `(-∞, +∞)`.
        public var range: (lower: T, upper: T) { dyn.range }

        // MARK: PDF/CDF/SF/Quantile

        /// The probability density function (PDF) evaluated at `x`.
        ///
        /// - Parameter x: Evaluation point (any real).
        /// - Returns: The density value `f(x)`.
        /// - Throws: Rethrows any underlying backend error.
        public func pdf(_ x: T) throws -> T { try dyn.pdf(x) }

        /// The natural logarithm of the PDF evaluated at `x`.
        ///
        /// - Parameter x: Evaluation point (any real).
        /// - Returns: The log-density `log f(x)`.
        /// - Throws: Rethrows any underlying backend error.
        public func logPdf(_ x: T) throws -> T { try dyn.logPdf(x) }

        /// The cumulative distribution function (CDF) evaluated at `x`.
        ///
        /// - Parameter x: Evaluation point (any real).
        /// - Returns: `P(X ≤ x)`.
        /// - Throws: Rethrows any underlying backend error.
        public func cdf(_ x: T) throws -> T { try dyn.cdf(x) }

        /// The survival function (SF) evaluated at `x`.
        ///
        /// - Parameter x: Evaluation point (any real).
        /// - Returns: `P(X > x)`.
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

        /// The mean of the Cauchy distribution.
        ///
        /// - Note: The mean does not exist (is undefined) for Cauchy. This property is therefore `nil`.
        public var mean: T? { nil }

        /// The variance of the distribution, if defined.
        ///
        /// - Note: For Cauchy the variance does not exist (is undefined). The backend reports `nil`.
        public var variance: T? { dyn.variance }

        /// The mode of the distribution, if defined.
        ///
        /// - Note: For Cauchy the mode exists and equals the location parameter.
        public var mode: T? { dyn.mode }

        /// The median of the distribution.
        ///
        /// - Note: For Cauchy the median equals the location parameter.
        public var median: T { dyn.median }

        /// The skewness of the distribution, if defined.
        ///
        /// - Note: For Cauchy the skewness is undefined; the backend reports `nil`.
        public var skewness: T? { dyn.skewness }

        /// The kurtosis of the distribution (Pearson), if defined.
        ///
        /// - Note: For Cauchy the kurtosis is undefined; the backend reports `nil`.
        public var kurtosis: T? { dyn.kurtosis }

        /// The excess kurtosis of the distribution, if defined.
        ///
        /// - Note: For Cauchy the excess kurtosis is undefined; the backend reports `nil`.
        public var kurtosisExcess: T? { dyn.kurtosisExcess }

        // MARK: Hazards

        /// The (instantaneous) hazard function `h(x) = f(x) / S(x)` evaluated at `x`.
        ///
        /// - Parameter x: Evaluation point (any real).
        /// - Returns: The hazard at `x`.
        /// - Throws: Rethrows any underlying backend error.
        public func hazard(_ x: T) throws -> T { try dyn.hazard(x) }

        /// The cumulative hazard function `H(x) = −log S(x)` evaluated at `x`.
        ///
        /// - Parameter x: Evaluation point (any real).
        /// - Returns: The cumulative hazard at `x`.
        /// - Throws: Rethrows any underlying backend error.
        public func chf(_ x: T) throws -> T { try dyn.chf(x) }

        // Lattice/discrete-only properties (continuous ⇒ nil)

        /// For lattice (discrete) distributions, the spacing between adjacent support points.
        ///
        /// - Note: Cauchy is continuous; this is `nil`.
        public var latticeStep: T? { nil }

        /// For lattice (discrete) distributions, the origin (offset) of the lattice.
        ///
        /// - Note: Cauchy is continuous; this is `nil`.
        public var latticeOrigin: T? { nil }

        /// The (differential) entropy h[X] where defined.
        ///
        /// - Note: For Cauchy with scale γ, h = ln(4πγ). The backend provides this value.
        public var entropy: T? { dyn.entropy }
    }
}
