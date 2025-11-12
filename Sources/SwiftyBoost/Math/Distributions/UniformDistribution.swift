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
    /// Continuous uniform (rectangular) distribution on `[lower, upper]`.
    public struct Uniform<T: Real & BinaryFloatingPoint & Sendable>: Sendable, DistributionProtocol {
        public typealias RealType = T

        /// Inclusive lower bound of the support interval.
        public let lower: T
        /// Inclusive upper bound of the support interval.
        public let upper: T

        /// Backing dynamic distribution that funnels calls into Boost.
        private let dyn: Distribution.Dynamic<T>

        /// Creates a continuous uniform distribution.
        ///
        /// - Parameters:
        ///   - lower: Lower bound (inclusive).
        ///   - upper: Upper bound (inclusive). Must satisfy `upper > lower`.
        public init(lower: T = 0, upper: T = 1) throws {
            guard upper > lower else {
                throw DistributionError.parameterOutOfRange(
                    name: "upper",
                    min: Double(lower),
                    max: Double.infinity
                )
            }
            self.lower = lower
            self.upper = upper
            self.dyn = try Distribution.Dynamic<T>(
                distributionName: "uniform",
                parameters: [
                    "lower": lower,
                    "upper": upper
                ]
            )
        }

        // MARK: Support

        public var supportLowerBound: T { dyn.supportLowerBound }
        public var supportUpperBound: T { dyn.supportUpperBound }
        public var range: (lower: T, upper: T) { dyn.range }

        // MARK: Core functions

        public func pdf(_ x: T) throws -> T { try dyn.pdf(x) }
        public func logPdf(_ x: T) throws -> T { try dyn.logPdf(x) }
        public func cdf(_ x: T) throws -> T { try dyn.cdf(x) }
        public func sf(_ x: T) throws -> T { try dyn.sf(x) }
        public func hazard(_ x: T) throws -> T { try dyn.hazard(x) }
        public func chf(_ x: T) throws -> T { try dyn.chf(x) }

        // MARK: Inverses

        public func quantile(_ p: T) throws -> T {
            guard p >= 0 && p <= 1 else {
                throw DistributionError.parameterOutOfRange(
                    name: "p",
                    min: 0.0,
                    max: 1.0
                )
            }
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

        public func klDivergence<D>(
            relativeTo other: D,
            options: Distribution.KLDivergenceOptions<T>
        ) throws -> T? where D: DistributionProtocol, D.RealType == T {
            try DistributionKLDivergenceHelper.evaluate(lhs: self, rhs: other, options: options)
        }
}
}
