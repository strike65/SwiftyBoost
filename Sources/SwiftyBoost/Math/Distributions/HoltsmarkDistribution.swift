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
    /// Holtsmark probability distribution (location–scale form).
    ///
    /// The Holtsmark distribution is a heavy–tailed, symmetric, continuous
    /// probability distribution that is strictly stable. This implementation
    /// exposes a location–scale parametrization:
    ///
    /// - `location` shifts the distribution along the real line.
    /// - `scale` stretches or shrinks the distribution and must be strictly positive.
    ///
    /// Core operations (PDF/CDF/log-PDF/SF/quantiles/hazards) delegate to a
    /// dynamic backend (`Distribution.Dynamic`) keyed by the distribution name
    /// `"holtsmark"`. This type is `Sendable` and suitable for use across
    /// concurrency domains.
    ///
    /// - Generic parameter `T`: The floating-point type used by this distribution
    ///   (`Float`, `Double`, or platform-available `Float80`). Must conform to
    ///   `Real & BinaryFloatingPoint & Sendable`.
    public struct Holtsmark<T: Real & BinaryFloatingPoint & Sendable>: Sendable,
        DistributionProtocol
    {
        /// The concrete real number type used by this distribution.
        typealias RealType = T

        /// Location parameter (shift).
        public let location: T

        /// Scale parameter (stretch). Must be strictly positive.
        public let scale: T

        /// Dynamic backend providing the numerical implementation for distribution functions.
        private let dyn: Distribution.Dynamic<T>

        /// Creates a Holtsmark distribution with the given parameters.
        ///
        /// - Parameters:
        ///   - loc: Location parameter. Defaults to `0`.
        ///   - scale: Scale parameter. Must be strictly positive. Defaults to `1`.
        ///
        /// - Throws: `DistributionError.parameterNotPositive(name: "scale", value: scale)`
        ///   if `scale <= 0`.
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
                distributionName: "holtsmark",
                parameters: [
                    "loc": loc,
                    "scale": scale
                ]
            )
        }

        // MARK: DistributionProtocol — Support

        /// Lower bound of the support.
        ///
        /// Delegates to the dynamic backend.
        public var supportLowerBound: T { dyn.supportLowerBound }

        /// Upper bound of the support.
        ///
        /// Delegates to the dynamic backend.
        public var supportUpperBound: T { dyn.supportUpperBound }

        /// Closed range of the support as a tuple `(lower, upper)`.
        ///
        /// Delegates to the dynamic backend.
        public var range: (lower: T, upper: T) { dyn.range }

        // MARK: PDF/CDF/SF/Quantile

        /// Probability density function evaluated at `x`.
        ///
        /// - Parameter x: The evaluation point.
        /// - Returns: The density value `f(x)`.
        /// - Throws: Rethrows any underlying backend error.
        public func pdf(_ x: T) throws -> T { try dyn.pdf(x) }

        /// Natural logarithm of the probability density function evaluated at `x`.
        ///
        /// - Parameter x: The evaluation point.
        /// - Returns: The log-density value `log f(x)`.
        /// - Throws: Rethrows any underlying backend error.
        public func logPdf(_ x: T) throws -> T { try dyn.logPdf(x) }

        /// Cumulative distribution function evaluated at `x`.
        ///
        /// - Parameter x: The evaluation point.
        /// - Returns: The CDF value `F(x)`.
        /// - Throws: Rethrows any underlying backend error.
        public func cdf(_ x: T) throws -> T { try dyn.cdf(x) }

        /// Survival function (complementary CDF) evaluated at `x`.
        ///
        /// - Parameter x: The evaluation point.
        /// - Returns: The SF value `S(x) = 1 - F(x)`.
        /// - Throws: Rethrows any underlying backend error.
        public func sf(_ x: T) throws -> T { try dyn.sf(x) }

        /// Quantile (inverse CDF) at probability `p`.
        ///
        /// - Parameter p: Probability in `[0, 1]`.
        /// - Returns: The value `x` such that `F(x) = p`.
        /// - Throws: `DistributionError` if `p` is out of bounds or the backend fails.
        public func quantile(_ p: T) throws -> T { try dyn.quantile(p) }

        /// Complement quantile (inverse survival function) at tail probability `q`.
        ///
        /// - Parameter q: Upper-tail probability in `[0, 1]`.
        /// - Returns: The value `x` such that `S(x) = q`.
        /// - Throws: `DistributionError` if `q` is out of bounds or the backend fails.
        public func quantileComplement(_ q: T) throws -> T {
            try dyn.quantileComplement(q)
        }

        // MARK: Moments

        /// Mean of the distribution, if it exists; otherwise `nil`.
        ///
        /// - Note: For heavy–tailed stable distributions, the mean can be undefined.
        ///   This implementation returns `0` for certain parameterizations by symmetry
        ///   and `nil` otherwise.
        public var mean: T? { dyn.mean }

        /// Variance of the distribution, if it exists; otherwise `nil`.
        ///
        /// Delegates to the dynamic backend.
        public var variance: T? { dyn.variance }

        /// Mode of the distribution, if defined; otherwise `nil`.
        ///
        /// Delegates to the dynamic backend.
        public var mode: T? { dyn.mode }

        /// Median of the distribution.
        ///
        /// Delegates to the dynamic backend.
        public var median: T { dyn.median }

        /// Skewness of the distribution, if it exists; otherwise `nil`.
        ///
        /// Delegates to the dynamic backend.
        public var skewness: T? { dyn.skewness }

        /// Kurtosis of the distribution (Pearson), if it exists; otherwise `nil`.
        ///
        /// Delegates to the dynamic backend.
        public var kurtosis: T? { dyn.kurtosis }

        /// Excess kurtosis of the distribution, if it exists; otherwise `nil`.
        ///
        /// Delegates to the dynamic backend.
        public var kurtosisExcess: T? { dyn.kurtosisExcess }

        // MARK: Hazards

        /// Hazard rate at `x`, if defined.
        ///
        /// - Parameter x: The evaluation point.
        /// - Returns: The hazard value `h(x) = f(x) / S(x)`.
        /// - Throws: Rethrows any underlying backend error.
        public func hazard(_ x: T) throws -> T { try dyn.hazard(x) }

        /// Cumulative hazard function at `x`, if defined.
        ///
        /// - Parameter x: The evaluation point.
        /// - Returns: The cumulative hazard `H(x) = -log S(x)`.
        /// - Throws: Rethrows any underlying backend error.
        public func chf(_ x: T) throws -> T { try dyn.chf(x) }

        // Lattice/discrete-only properties (continuous ⇒ nil)

        /// Step size for lattice (discrete) distributions.
        ///
        /// Holtsmark is continuous, so this is `nil`.
        public var latticeStep: T? { nil }

        /// Origin for lattice (discrete) distributions.
        ///
        /// Holtsmark is continuous, so this is `nil`.
        public var latticeOrigin: T? { nil }

        /// Differential entropy, if defined.
        ///
        /// For location–scale families, the entropy decomposes as
        /// `H(X) = log(scale) + H(Z)` where `Z` is the standardized distribution
        /// with `scale = 1` and `location = 0`. The constant term for the
        /// standardized Holtsmark distribution is provided by
        /// `Constants.holtsmarkEntropy()`.
        ///
        /// - Returns: The entropy in nats as `T`, or `nil` if it cannot be computed.
        public var entropy: T? {
            return T.log(self.scale) + Constants.holtsmarkEntropy()
        }
    }
}
