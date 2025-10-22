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
import CBoostBridge
import Foundation

extension Distribution {
    /// A continuous arcsine distribution on a finite interval.
    ///
    /// This distribution is supported on the closed interval `[minX, maxX]` and has
    /// the probability density function (PDF)
    ///
    ///     f(x) = 1 / (π * sqrt((x - minX) * (maxX - x)))   for x ∈ (minX, maxX)
    ///
    /// The distribution is U-shaped, with singularities at the endpoints of the support.
    ///
    /// - Note: This type is generic over a floating-point type `T` and supports `Float`, `Double`,
    ///   and, on x86 architectures, `Float80`. Internally it bridges to Boost via `CBoostBridge`.
    public struct Arcsine<T: BinaryFloatingPoint & Sendable>: Sendable, DistributionProtocol {
        /// The floating-point type used by this distribution.
        typealias Real = T

        /// The lower bound of the support (a).
        public let minX: T
        /// The upper bound of the support (b).
        public let maxX: T

        private let dyn: Distribution.Dynamic<T>

        /// Creates an arcsine distribution with the given support.
        ///
        /// - Parameters:
        ///   - min_x: The lower bound `a` of the support. Must be finite and strictly less than `max_x`.
        ///   - max_x: The upper bound `b` of the support. Must be finite and strictly greater than `min_x`.
        ///
        /// - Throws: `DistributionError.invalidCombination` if `min_x >= max_x`,
        ///           `DistributionError.parameterNotFinite` if either bound is not finite,
        ///           or `DistributionError.generalError` if the underlying distribution could not be initialized.
        public init(minX min_x: T = 0, maxX max_x: T = 0) throws {
            guard min_x < max_x else {
                throw DistributionError.invalidCombination(message: "minX must be less than maxX")
            }
            guard min_x.isFinite else {
                throw DistributionError.parameterNotFinite(name: "minX")
            }
            guard max_x.isFinite else {
                throw DistributionError.parameterNotFinite(name: "maxX")
            }

            self.minX = min_x
            self.maxX = max_x
            self.dyn = try Distribution.Dynamic<T>(
                distributionName: "arcsine",
                parameters: [
                    "minX": min_x,
                    "maxX": max_x,
                ]
            )
        }

        /// The lower bound of the support.
        public var supportLowerBound: T { dyn.supportLowerBound }

        /// The upper bound of the support.
        public var supportUpperBound: T { dyn.supportUpperBound }

        /// The support as a (lower, upper) tuple.
        public var range: (lower: T, upper: T) { dyn.range }

        /// The probability density function (PDF).
        ///
        /// - Parameter x: Evaluation point.
        /// - Returns: The density at `x`. Values outside the support yield zero.
        public func pdf(_ x: T) throws -> T { try dyn.pdf(x) }

        /// The natural logarithm of the PDF.
        ///
        /// - Parameter x: Evaluation point.
        /// - Returns: `log(pdf(x))`. May be `-infinity` at the boundaries.
        public func logPdf(_ x: T) throws -> T { try dyn.logPdf(x) }

        /// The cumulative distribution function (CDF).
        ///
        /// - Parameter x: Evaluation point.
        /// - Returns: `P(X ≤ x)`.
        public func cdf(_ x: T) throws -> T { try dyn.cdf(x) }

        /// The survival function (SF), i.e. the complementary CDF.
        ///
        /// - Parameter x: Evaluation point.
        /// - Returns: `P(X > x) = 1 - CDF(x)`.
        public func sf(_ x: T) throws -> T { try dyn.sf(x) }

        /// The quantile function (inverse CDF).
        ///
        /// - Parameter p: A probability in `[0, 1]`.
        /// - Returns: `x` such that `P(X ≤ x) = p`.
        public func quantile(_ p: T) throws -> T { try dyn.quantile(p) }

        /// The complementary quantile function (inverse survival function).
        ///
        /// - Parameter q: A probability in `[0, 1]`.
        /// - Returns: `x` such that `P(X > x) = q`.
        public func quantileComplement(_ q: T) throws -> T { try dyn.quantileComplement(q) }

        /// The mean of the distribution.
        ///
        /// - Note: Returns `nil` if not finite for the chosen numeric type.
        public var mean: T? { dyn.mean }

        /// The variance of the distribution.
        ///
        /// - Note: Returns `nil` if not finite for the chosen numeric type.
        public var variance: T? { dyn.variance }

        /// The mode of the distribution.
        ///
        /// - Note: For the arcsine distribution on a finite interval, the density is unbounded at the endpoints.
        ///   This value reflects the underlying implementation choice in Boost.
        public var mode: T? { dyn.mode }

        /// The median of the distribution.
        public var median: T { dyn.median }

        /// The skewness of the distribution.
        ///
        /// - Note: Returns `nil` if not finite for the chosen numeric type.
        public var skewness: T? { dyn.skewness }

        /// The kurtosis of the distribution.
        ///
        /// - Note: Returns `nil` if not finite for the chosen numeric type.
        public var kurtosis: T? { dyn.kurtosis }

        /// The excess kurtosis of the distribution (`kurtosis - 3`).
        ///
        /// - Note: Returns `nil` if not finite for the chosen numeric type.
        public var kurtosisExcess: T? { dyn.kurtosisExcess }

        /// The hazard (failure) rate function `h(x) = f(x) / (1 - F(x))`.
        ///
        /// - Parameter x: Evaluation point.
        /// - Returns: The hazard at `x`.
        public func hazard(_ x: T) throws -> T { try dyn.hazard(x) }

        /// The cumulative hazard function `H(x) = -log(1 - F(x))`.
        ///
        /// - Parameter x: Evaluation point.
        /// - Returns: The cumulative hazard at `x`.
        public func chf(_ x: T) throws -> T { try dyn.chf(x) }

        /// The lattice step size, if the distribution is lattice-supported.
        ///
        /// Always `nil` for this continuous distribution.
        public var latticeStep: T? { nil }

        /// The lattice origin, if the distribution is lattice-supported.
        ///
        /// Always `nil` for this continuous distribution.
        public var latticeOrigin: T? { nil }

        /// The differential entropy of the arcsine distribution.
        ///
        /// For support width `w = maxX - minX`, the entropy is:
        ///
        ///     H = ln(π / 4) + ln(w) = ln(π w / 4)
        ///
        /// - Returns: The entropy as `T`, or `nil` if not finite for the chosen numeric type.
        public var entropy: T? {
            guard self.minX < self.maxX else { return nil }
            let width = self.maxX - self.minX
            #if arch(x86_64) || arch(i386)
            let qpi: Float80 = Constants.quarterPi
            return T(log(qpi)) + T(log(Float80(width)))
            #else
            let qpi: Double = Constants.quarterPi
            return T(log(qpi)) + T(log(Double(width)))
            #endif
        }
    }
}
