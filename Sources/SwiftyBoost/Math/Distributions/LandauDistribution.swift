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
    /// Landau distribution (location–scale parameterisation).
    ///
    /// The Landau distribution is a heavy-tailed, strictly stable distribution
    /// that arises in modelling energy loss of charged particles traversing thin
    /// materials. Only the location and scale parameters are finite; all algebraic
    /// moments diverge.
    public struct Landau<T: Real & BinaryFloatingPoint & Sendable>: Sendable, DistributionProtocol {
        public typealias RealType = T

        /// Location (shift) parameter.
        public let location: T

        /// Scale parameter (must be strictly positive).
        public let scale: T

        /// Dynamic backend used for numerical evaluations.
        private let dyn: Distribution.Dynamic<T>

        /// Creates a Landau distribution with the specified parameters.
        ///
        /// - Parameters:
        ///   - location: Location (shift) parameter. Defaults to `0`.
        ///   - scale: Scale parameter. Must be strictly positive. Defaults to `1`.
        /// - Throws: ``DistributionError/parameterNotPositive`` when `scale <= 0`.
        public init(location: T = 0, scale: T = 1) throws {
            guard scale > 0 else {
                throw DistributionError.parameterNotPositive(
                    name: "scale",
                    value: scale
                )
            }
            self.location = location
            self.scale = scale
            self.dyn = try Distribution.Dynamic<T>(
                distributionName: "landau",
                parameters: [
                    "location": location,
                    "scale": scale
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

        // MARK: - Summary statistics

        /// The Landau distribution has no finite mean.
        public var mean: T? { nil }

        /// The Landau distribution has no finite variance.
        public var variance: T? { nil }

        public var mode: T? { dyn.mode }
        public var median: T { dyn.median }

        /// Higher moments diverge and are therefore reported as `nil`.
        public var skewness: T? { nil }
        public var kurtosis: T? { nil }
        public var kurtosisExcess: T? { nil }

        public var entropy: T? { dyn.entropy }

        // MARK: - Lattice metadata

        public var latticeStep: T? { nil }
        public var latticeOrigin: T? { nil }
        public var isDiscrete: Bool { dyn.isDiscrete }

        /// Estimates the Kullback–Leibler divergence numerically using the dynamic backend.
        public func klDivergence<D>(
            relativeTo other: D,
            options: Distribution.KLDivergenceOptions<T>
        ) throws -> T? where D: DistributionProtocol, D.RealType == T {
            try DistributionKLDivergenceHelper.evaluate(lhs: self, rhs: other, options: options)
        }
}
}
