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
    /// Non-central Fisher–Snedecor F distribution backed by Boost.Math.
    public struct NonCentralF<T: Real & BinaryFloatingPoint & Sendable>: Sendable, DistributionProtocol {

        public typealias RealType = T

        /// Numerator degrees of freedom (ν₁ > 0).
        public let degreesOfFreedom1: T

        /// Denominator degrees of freedom (ν₂ > 0).
        public let degreesOfFreedom2: T

        /// Non-centrality parameter (λ ≥ 0).
        public let nonCentrality: T

        private let dyn: Distribution.Dynamic<T>

        /// Creates a non-central F distribution `F(ν₁, ν₂, λ)`.
        public init(degreesOfFreedom1 v1: T,
                    degreesOfFreedom2 v2: T,
                    nonCentrality lambda: T) throws {
            guard v1.isFinite else {
                throw DistributionError<T>.parameterNotFinite(name: "degreesOfFreedom1", value: v1)
            }
            guard v2.isFinite else {
                throw DistributionError<T>.parameterNotFinite(name: "degreesOfFreedom2", value: v2)
            }
            guard lambda.isFinite else {
                throw DistributionError<T>.parameterNotFinite(name: "nonCentrality", value: lambda)
            }
            guard v1 > 0 else {
                throw DistributionError<T>.parameterNotPositive(name: "degreesOfFreedom1", value: v1)
            }
            guard v2 > 0 else {
                throw DistributionError<T>.parameterNotPositive(name: "degreesOfFreedom2", value: v2)
            }
            guard lambda >= 0 else {
                throw DistributionError<T>.parameterOutOfRange(name: "nonCentrality", min: 0, max: .infinity)
            }

            self.degreesOfFreedom1 = v1
            self.degreesOfFreedom2 = v2
            self.nonCentrality = lambda
            self.dyn = try Distribution.Dynamic<T>(
                distributionName: "non_central_f",
                parameters: [
                    "df1": v1,
                    "df2": v2,
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

        /// The mean `E[X]`, available when `ν₂ > 2`.
        public var mean: T? {
            if degreesOfFreedom2 > 2 {
                let numerator = degreesOfFreedom2 * (degreesOfFreedom1 + nonCentrality)
                let denominator = degreesOfFreedom1 * (degreesOfFreedom2 - 2)
                return numerator / denominator
            }
            return nil
        }

        /// The variance `Var[X]`, available when `ν₂ > 4`.
        public var variance: T? {
            if degreesOfFreedom2 > 4 {
                let v1 = degreesOfFreedom1
                let v2 = degreesOfFreedom2
                let lambda = nonCentrality
                let ratio = v2 / v1
                let ratioSquared = ratio * ratio
                let term1 = (v1 + lambda) * (v1 + lambda)
                let term2 = (v1 + T(2) * lambda) * (v2 - 2)
                let numerator = T(2) * ratioSquared * (term1 + term2)
                let denominator = (v2 - 2) * (v2 - 2) * (v2 - 4)
                return numerator / denominator
            }
            return nil
        }

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

        public func klDivergence(
            relativeTo other: Self,
            options: Distribution.KLDivergenceOptions<T> = .automatic()
        ) throws -> T? {
            try dyn.klDivergence(relativeTo: other.dyn, options: options)
        }
    }
}
