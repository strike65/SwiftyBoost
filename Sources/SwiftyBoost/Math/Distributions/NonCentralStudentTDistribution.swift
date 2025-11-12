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
    /// Non-central Student’s t distribution backed by Boost.Math.
    public struct NonCentralStudentT<T: Real & BinaryFloatingPoint & Sendable>: Sendable, DistributionProtocol {

        public typealias RealType = T

        /// Degrees of freedom ν (> 0).
        public let degreesOfFreedom: T

        /// Non-centrality parameter λ (finite real).
        public let nonCentrality: T

        private let dyn: Distribution.Dynamic<T>

        /// Creates a non-central Student’s t distribution `t(ν, λ)`.
        public init(degreesOfFreedom v: T, nonCentrality lambda: T) throws {
            guard v.isFinite else {
                throw DistributionError<T>.parameterNotFinite(name: "degreesOfFreedom", value: v)
            }
            guard lambda.isFinite else {
                throw DistributionError<T>.parameterNotFinite(name: "nonCentrality", value: lambda)
            }
            guard v > 0 else {
                throw DistributionError<T>.parameterNotPositive(name: "degreesOfFreedom", value: v)
            }

            self.degreesOfFreedom = v
            self.nonCentrality = lambda
            self.dyn = try Distribution.Dynamic<T>(
                distributionName: "non_central_t",
                parameters: [
                    "df": v,
                    "lambda": lambda
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
        public func quantile(_ p: T) throws -> T { try dyn.quantile(p) }
        public func quantileComplement(_ q: T) throws -> T { try dyn.quantileComplement(q) }

        // MARK: - Moments

        public var mean: T? { dyn.mean }
        public var variance: T? { dyn.variance }
        public var mode: T? { dyn.mode }
        public var median: T { dyn.median }
        public var skewness: T? { dyn.skewness }
        public var kurtosis: T? { dyn.kurtosis }
        public var kurtosisExcess: T? { dyn.kurtosisExcess }

        // MARK: - Hazards

        public func hazard(_ x: T) throws -> T { try dyn.hazard(x) }
        public func chf(_ x: T) throws -> T { try dyn.chf(x) }

        // MARK: - Lattice metadata

        public var latticeStep: T? { nil }
        public var latticeOrigin: T? { nil }

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
