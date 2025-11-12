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
    /// Poisson distribution backed by Boost.Math.
    public struct Poisson<T: Real & BinaryFloatingPoint & Sendable>: Sendable, DistributionProtocol {

        public typealias RealType = T

        /// Mean `λ ≥ 0`.
        public let lambda: T

        private let dyn: Distribution.Dynamic<T>

        /// Creates a Poisson distribution with the given mean.
        public init(lambda: T) throws {
            guard lambda.isFinite else {
                throw DistributionError<T>.parameterNotFinite(name: "lambda", value: lambda)
            }
            guard lambda >= 0 else {
                throw DistributionError<T>.parameterOutOfRange(name: "lambda", min: 0, max: .infinity)
            }
            self.lambda = lambda
            self.dyn = try Distribution.Dynamic<T>(
                distributionName: "poisson",
                parameters: ["lambda": lambda]
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
        public func quantile(_ p: T) throws -> T { try dyn.quantile(p) }
        public func quantileComplement(_ q: T) throws -> T { try dyn.quantileComplement(q) }

        // MARK: - Moments

        public var mean: T? { lambda }

        public var variance: T? { lambda }

        /// Mode is ⌊λ⌋ (Boost returns the principal mode).
        public var mode: T? { dyn.mode }

        public var median: T { dyn.median }

        /// Skewness = 1 / √λ (undefined at λ = 0).
        public var skewness: T? {
            if lambda > 0 {
                return 1 / lambda.squareRoot()
            }
            return nil
        }

        /// Kurtosis (Pearson’s β₂) = 3 + 1/λ.
        public var kurtosis: T? {
            if lambda > 0 {
                return 3 + 1 / lambda
            }
            return nil
        }

        /// Excess kurtosis = 1/λ.
        public var kurtosisExcess: T? {
            if lambda > 0 {
                return 1 / lambda
            }
            return nil
        }

        // MARK: - Hazards

        public func hazard(_ x: T) throws -> T { try dyn.hazard(x) }
        public func chf(_ x: T) throws -> T { try dyn.chf(x) }

        // MARK: - Lattice metadata

        public var latticeStep: T? { 1 }
        public var latticeOrigin: T? { 0 }

        // MARK: - Entropy & divergence

        public var entropy: T? { dyn.entropy }
        public var isDiscrete: Bool { dyn.isDiscrete }

        public func klDivergence<D>(
            relativeTo other: D,
            options: Distribution.KLDivergenceOptions<T>
        ) throws -> T? where D: DistributionProtocol, D.RealType == T {
            try DistributionKLDivergenceHelper.evaluate(lhs: self, rhs: other, options: options)
        }
}
}
