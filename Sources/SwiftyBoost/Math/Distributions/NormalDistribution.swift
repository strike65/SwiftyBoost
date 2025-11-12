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
    public struct Normal<T: Real & BinaryFloatingPoint & Sendable>: Sendable, DistributionProtocol {
        public typealias RealType = T
        public let location: T
        public let sd: T
        private let dyn: Distribution.Dynamic<T>
        public init(mean loc: T = 0, sd: T = 1) throws {
            guard sd > 0 else {
                throw DistributionError.parameterNotPositive(
                    name: "sd",
                    value: sd
                )
            }
            self.location = loc
            self.sd = sd
            self.dyn = try Distribution.Dynamic<T>(
                distributionName: "normal",
                parameters: [
                    "mean": loc,
                    "sd": sd
                ]
            )
        }
        /// Lower bound of the support (`−∞` in the ideal Gaussian).
        public var supportLowerBound: T { dyn.supportLowerBound }
        /// Upper bound of the support (`+∞` in the ideal Gaussian).
        public var supportUpperBound: T { dyn.supportUpperBound }
        /// Convenience tuple containing both support bounds.
        public var range: (lower: T, upper: T) { dyn.range }
        /// Probability density function evaluated at `x`.
        public func pdf(_ x: T) throws -> T { try dyn.pdf(x) }
        /// Natural logarithm of the PDF evaluated at `x`.
        public func logPdf(_ x: T) throws -> T { try dyn.logPdf(x) }
        /// Cumulative distribution function `F(x) = P(X ≤ x)`.
        public func cdf(_ x: T) throws -> T { try dyn.cdf(x) }
        /// Survival function `S(x) = P(X > x)`.
        public func sf(_ x: T) throws -> T { try dyn.sf(x) }
        /// Lower-tail quantile (inverse CDF).
        public func quantile(_ p: T) throws -> T { try dyn.quantile(p) }
        /// Upper-tail quantile (inverse survival function).
        public func quantileComplement(_ q: T) throws -> T {
            try dyn.quantileComplement(q)
        }
        public var mean: T? { dyn.mean }
        public var variance: T? { dyn.variance }
        public var mode: T? { dyn.mode }
        public var median: T { dyn.median }
        public var skewness: T? { dyn.skewness }
        public var kurtosis: T? { dyn.kurtosis }
        public var kurtosisExcess: T? { dyn.kurtosisExcess }
        /// Instantaneous hazard function `h(x) = f(x) / S(x)`.
        public func hazard(_ x: T) throws -> T { try dyn.hazard(x) }
        /// Cumulative hazard function `H(x) = −log S(x)`.
        public func chf(_ x: T) throws -> T { try dyn.chf(x) }
        /// Lattice spacing for discrete distributions (not applicable to Gaussian).
        public var latticeStep: T? { nil }
        /// Lattice origin for discrete distributions (not applicable to Gaussian).
        public var latticeOrigin: T? { nil }
        /// Differential entropy `h[X] = ½ ln(2πeσ²)` supplied by the backend.
        public var entropy: T? { dyn.entropy }

        /// Indicates whether this distribution is discrete (`true`) or continuous (`false`).
        ///
        /// Normal (Gaussian) distributions are continuous, so this always returns `false`.
        public var isDiscrete: Bool { dyn.isDiscrete }

        /// Computes the Kullback–Leibler divergence `D_KL(self || other)` using the analytic closed form.
        ///
        /// - Parameters:
        ///   - other: The reference normal distribution *Q*.
        ///   - options: Unused; included for signature compatibility. You can pass ``Distribution/KLDivergenceOptions/automatic()``.
        /// - Returns: The divergence in nats, or `nil` if either standard deviation is non-positive.
        /// - Throws: Never throws directly; kept for parity with other distributions.
        public func klDivergence<D: DistributionProtocol>(
            relativeTo other: D,
            options: Distribution.KLDivergenceOptions<T>
        ) throws -> T? where D.RealType == T {
            if let rhs = other as? Self {
                let sigmaSelf = self.sd
                let sigmaOther = rhs.sd
                guard sigmaSelf > 0, sigmaOther > 0 else { return nil }
                let logTerm = T.log(sigmaOther / sigmaSelf)
                let meanDiff = self.location - rhs.location
                let varianceTerm = (sigmaSelf * sigmaSelf + meanDiff * meanDiff) / (2 * sigmaOther * sigmaOther)
                return logTerm + varianceTerm - 0.5
            }
            return try DistributionKLDivergenceHelper.evaluate(lhs: self, rhs: other, options: options)
        }
    }
}
