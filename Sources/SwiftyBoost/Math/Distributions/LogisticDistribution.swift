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
    /// Logistic distribution (location-scale form).
    ///
    /// The logistic distribution is a symmetric, heavy-tailed alternative to the normal
    /// distribution. It is parameterised by a location parameter `μ` and a strictly
    /// positive scale parameter `s`.
    public struct Logistic<T: Real & BinaryFloatingPoint & Sendable>: Sendable, DistributionProtocol {
        typealias RealType = T

        /// Location parameter.
        public let location: T

        /// Scale parameter (must be strictly positive).
        public let scale: T

        private let dyn: Distribution.Dynamic<T>

        /// Creates a logistic distribution with the provided parameters.
        ///
        /// - Parameters:
        ///   - location: Location parameter `μ`. Defaults to `0`.
        ///   - scale: Scale parameter `s`. Must be strictly positive. Defaults to `1`.
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
                distributionName: "logistic",
                parameters: [
                    "location": location,
                    "scale": scale
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

        public func quantile(_ p: T) throws -> T { try dyn.quantile(p) }
        public func quantileComplement(_ q: T) throws -> T { try dyn.quantileComplement(q) }

        // MARK: Summary statistics

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

        /// Computes the Kullback–Leibler divergence numerically via the dynamic backend.
        public func klDivergence(
            relativeTo other: Self,
            options: Distribution.KLDivergenceOptions<T>
        ) throws -> T? {
            try dyn.klDivergence(relativeTo: other.dyn, options: options)
        }
    }
}
