//
//  Created by Volker Thieme
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
    /// Rayleigh distribution with scale parameter σ > 0.
    ///
    /// Definitions:
    /// - Support: `[0, +∞)`
    /// - PDF: `f(x) = (x / σ²) · exp(−x² / (2 σ²))` for `x ≥ 0`
    /// - CDF: `F(x) = 1 − exp(−x² / (2 σ²))`
    /// - Survival: `S(x) = exp(−x² / (2 σ²))`
    ///
    /// Moments:
    /// - Mean: `σ √(π / 2)`
    /// - Variance: `(2 − π / 2) σ²`
    /// - Mode: `σ`
    /// - Median: `σ √(2 ln 2)`
    /// - Skewness: `(2 √π (π − 3)) / (4 − π)^(3/2)`
    /// - Excess kurtosis: `(−6 π² + 24 π − 16) / (4 − π)²`
    /// - Entropy: `1 + ln(σ / √2) + γ / 2` (γ Euler–Mascheroni)
    public struct Rayleigh<T: Real & BinaryFloatingPoint & Sendable>: Sendable, DistributionProtocol {
        public typealias RealType = T

        /// Scale parameter σ (> 0).
        public let scale: T

        /// Backing dynamic distribution that performs all evaluations.
        private let dyn: Distribution.Dynamic<T>

        /// Creates a Rayleigh distribution with the provided scale parameter.
        ///
        /// - Parameter scale: σ, must be strictly positive.
        public init(scale: T) throws {
            guard scale > 0 else {
                throw DistributionError.parameterNotPositive(
                    name: "scale",
                    value: scale
                )
            }
            self.scale = scale
            self.dyn = try Distribution.Dynamic<T>(
                distributionName: "rayleigh",
                parameters: ["scale": scale]
            )
        }

        // MARK: Support

        public var supportLowerBound: T { dyn.supportLowerBound }
        public var supportUpperBound: T { dyn.supportUpperBound }
        public var range: (lower: T, upper: T) { dyn.range }

        // MARK: Core functions

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

        public func logPdf(_ x: T) throws -> T {
            guard x >= 0 else {
                throw DistributionError.parameterOutOfRange(
                    name: "x",
                    min: 0.0,
                    max: Double.infinity
                )
            }
            return try dyn.logPdf(x)
        }

        public func cdf(_ x: T) throws -> T {
            guard x >= 0 else {
                throw DistributionError.parameterOutOfRange(
                    name: "x",
                    min: 0.0,
                    max: Double.infinity
                )
            }
            return try dyn.cdf(x)
        }

        public func sf(_ x: T) throws -> T {
            guard x >= 0 else {
                throw DistributionError.parameterOutOfRange(
                    name: "x",
                    min: 0.0,
                    max: Double.infinity
                )
            }
            return try dyn.sf(x)
        }

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

        // MARK: Inverses

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

        public var mean: T? { dyn.mean }
        public var variance: T? { dyn.variance }
        public var mode: T? { dyn.mode }
        public var median: T { dyn.median }
        public var skewness: T? { dyn.skewness }
        public var kurtosis: T? { dyn.kurtosis }
        public var kurtosisExcess: T? { dyn.kurtosisExcess }
        public var entropy: T? { dyn.entropy }

        // MARK: Lattice metadata

        public var latticeStep: T? { nil }
        public var latticeOrigin: T? { nil }
        public var isDiscrete: Bool { dyn.isDiscrete }

        public func klDivergence(
            relativeTo other: Self,
            options: Distribution.KLDivergenceOptions<T>
        ) throws -> T? {
            try dyn.klDivergence(relativeTo: other.dyn, options: options)
        }
    }
}
