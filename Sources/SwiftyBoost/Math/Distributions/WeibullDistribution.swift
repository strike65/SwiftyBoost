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
    /// Weibull distribution with shape `k > 0` and scale `λ > 0`.
    ///
    /// Support: `[0, +∞)`. Includes the exponential distribution (k = 1) and
    /// approximations to Rayleigh (k = 2) and extreme value forms.
    public struct Weibull<T: Real & BinaryFloatingPoint & Sendable>: Sendable, DistributionProtocol {
        public typealias RealType = T

        /// Shape parameter `k` controlling the hazard profile.
        public let shape: T
        /// Scale parameter `λ` that stretches the distribution.
        public let scale: T

        /// Backing dynamic distribution that delegates to Boost.
        private let dyn: Distribution.Dynamic<T>

        /// Creates a Weibull distribution.
        ///
        /// - Parameters:
        ///   - shape: Shape parameter `k`, must be strictly positive.
        ///   - scale: Scale parameter `λ`, must be strictly positive (default 1).
        public init(shape: T, scale: T = 1) throws {
            guard shape > 0 else {
                throw DistributionError.parameterNotPositive(
                    name: "shape",
                    value: shape
                )
            }
            guard scale > 0 else {
                throw DistributionError.parameterNotPositive(
                    name: "scale",
                    value: scale
                )
            }
            self.shape = shape
            self.scale = scale
            self.dyn = try Distribution.Dynamic<T>(
                distributionName: "weibull",
                parameters: [
                    "shape": shape,
                    "scale": scale
                ]
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
