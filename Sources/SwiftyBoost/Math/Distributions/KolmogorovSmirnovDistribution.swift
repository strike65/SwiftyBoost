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
    /// Kolmogorov–Smirnov distribution for the one-sample test statistic.
    ///
    /// The Kolmogorov–Smirnov distribution models the limiting distribution of the
    /// supremum of the absolute difference between an empirical CDF and a reference
    /// CDF. This implementation follows the Boost.Math parametrisation with a single
    /// positive parameter `n` corresponding to the finite sample size.
    ///
    /// - Generic parameter `T`: Any floating-point type supported by ``Distribution/Dynamic``
    ///   (`Float`, `Double`, or `Float80` on supported platforms).
    public struct KolmogorovSmirnov<T: Real & BinaryFloatingPoint & Sendable>: Sendable, DistributionProtocol {
        public typealias RealType = T

        /// Number of observations in the empirical sample.
        public let numberOfObservations: T

        /// Dynamic backend that provides the numerical implementation.
        private let dyn: Distribution.Dynamic<T>

        /// Creates a Kolmogorov–Smirnov distribution for a sample of size `n`.
        ///
        /// - Parameter n: The (strictly positive) number of observations.
        /// - Throws: ``DistributionError/parameterNotPositive(name:value:)`` when `n <= 0`.
        public init(numberOfObservations n: T) throws {
            guard n > 0 else {
                throw DistributionError.parameterNotPositive(
                    name: "numberOfObservations",
                    value: n
                )
            }
            self.numberOfObservations = n
            self.dyn = try Distribution.Dynamic<T>(
                distributionName: "kolmogorov_smirnov",
                parameters: [
                    "n": n
                ]
            )
        }

        // MARK: - Support

        public var supportLowerBound: T { dyn.supportLowerBound }
        public var supportUpperBound: T { dyn.supportUpperBound }
        public var range: (lower: T, upper: T) { dyn.range }

        // MARK: - Core functions

        public func pdf(_ x: T) throws -> T { try dyn.pdf(x) }
        public func logPdf(_ x: T) throws -> T { try dyn.logPdf(x) }
        public func cdf(_ x: T) throws -> T { try dyn.cdf(x) }
        public func sf(_ x: T) throws -> T { try dyn.sf(x) }
        public func hazard(_ x: T) throws -> T { try dyn.hazard(x) }
        public func chf(_ x: T) throws -> T { try dyn.chf(x) }

        // MARK: - Inverses

        public func quantile(_ p: T) throws -> T { try dyn.quantile(p) }
        public func quantileComplement(_ q: T) throws -> T { try dyn.quantileComplement(q) }

        // MARK: - Moments and summary statistics

        public var mean: T? { dyn.mean }
        public var variance: T? { dyn.variance }
        public var mode: T? { dyn.mode }
        public var median: T { dyn.median }
        public var skewness: T? { dyn.skewness }
        public var kurtosis: T? { dyn.kurtosis }
        public var kurtosisExcess: T? { dyn.kurtosisExcess }
        public var entropy: T? { dyn.entropy }

        // MARK: - Lattice metadata

        public var latticeStep: T? { nil }
        public var latticeOrigin: T? { nil }
        public var isDiscrete: Bool { dyn.isDiscrete }

        /// Computes the Kullback–Leibler divergence against another Kolmogorov–Smirnov distribution.
        public func klDivergence<D>(
            relativeTo other: D,
            options: Distribution.KLDivergenceOptions<T>
        ) throws -> T? where D: DistributionProtocol, D.RealType == T {
            try DistributionKLDivergenceHelper.evaluate(lhs: self, rhs: other, options: options)
        }
}
}
