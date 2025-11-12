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

    /// Non-central chi-squared probability distribution.
    ///
    /// This distribution generalises the central Chi-squared law to the case where the
    /// squared Gaussian components have non-zero means. It is parameterised by the
    /// degrees of freedom `ν > 0` and the non-centrality parameter `λ ≥ 0`.
    ///
    /// - Generic parameter `T`: Any binary floating-point type supported by the runtime
    ///   (Double, Float, Float80 on supported architectures).
    ///
    /// - SeeAlso:
    ///   - Wikipedia: Noncentral chi-squared distribution
    ///   - Boost.Math `non_central_chi_squared_distribution`
    public struct NonCentralChiSquared<T: Real & BinaryFloatingPoint & Sendable>: Sendable, DistributionProtocol {

        /// Floating-point type used throughout the distribution.
        public typealias RealType = T

        /// Degrees of freedom `ν`.
        public let degreesOfFreedom: T

        /// Non-centrality parameter `λ`.
        public let nonCentrality: T

        /// Boost-backed dynamic delegate providing core functionality.
        private let dyn: Distribution.Dynamic<T>

        /// Creates a non-central chi-squared distribution.
        ///
        /// - Parameters:
        ///   - degreesOfFreedom: Positive degrees of freedom `ν`.
        ///   - nonCentrality: Non-centrality parameter `λ ≥ 0`.
        ///
        /// - Throws:
        ///   - ``DistributionError/parameterNotPositive`` if `degreesOfFreedom ≤ 0`.
        ///   - ``DistributionError/parameterOutOfRange`` if `nonCentrality < 0`.
        ///   - ``DistributionError/parameterNotFinite`` if either argument is non-finite.
        public init(degreesOfFreedom: T, nonCentrality lambda: T) throws {
            guard degreesOfFreedom.isFinite else {
                throw DistributionError<T>.parameterNotFinite(
                    name: "degreesOfFreedom",
                    value: degreesOfFreedom
                )
            }
            guard lambda.isFinite else {
                throw DistributionError<T>.parameterNotFinite(
                    name: "nonCentrality",
                    value: lambda
                )
            }
            guard degreesOfFreedom > 0 else {
                throw DistributionError<T>.parameterNotPositive(
                    name: "degreesOfFreedom",
                    value: degreesOfFreedom
                )
            }
            guard lambda >= 0 else {
                throw DistributionError<T>.parameterOutOfRange(
                    name: "nonCentrality",
                    min: 0,
                    max: .infinity
                )
            }

            self.degreesOfFreedom = degreesOfFreedom
            self.nonCentrality = lambda
            self.dyn = try Distribution.Dynamic<T>(
                distributionName: "non_central_chi_squared",
                parameters: [
                    "df": degreesOfFreedom,
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

        /// Mean `E[X] = ν + λ`.
        public var mean: T? { self.moments.mu1 }

        /// Variance `Var[X] = 2(ν + 2λ)`.
        public var variance: T? {
            let value = self.moments.mu2
            return value.isFinite ? value : nil
        }

        /// Mode supplied by the backend (no simple closed form).
        public var mode: T? { dyn.mode }

        /// Median supplied by the backend (numerical).
        public var median: T { dyn.median }

        /// Skewness computed from the third central moment.
        public var skewness: T? {
            let m = self.moments
            guard m.mu2 > 0 else { return nil }
            let sigma = m.mu2.squareRoot()
            let denom = m.mu2 * sigma
            guard denom > 0 else { return nil }
            return m.mu3 / denom
        }

        /// Kurtosis (Pearson’s β₂) computed from the fourth central moment.
        public var kurtosis: T? {
            let m = self.moments
            guard m.mu2 > 0 else { return nil }
            let denom = m.mu2 * m.mu2
            guard denom > 0 else { return nil }
            return m.mu4 / denom
        }

        /// Excess kurtosis `β₂ − 3`.
        public var kurtosisExcess: T? {
            if let kurt = self.kurtosis {
                return kurt - 3
            }
            return nil
        }

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
// helpers
        public static func find_degreesOfFreedom(
            lambda: T,
            x: T,
            p: T) -> Int {
            let res: T
            if T.self == Double.self {
                res = T(
                    bs_non_central_chisquare_find_degreesOfFreedom_d(
                        Double(lambda),
                        Double(x),
                        Double(p))
                    )
            }
            else if T.self == Float.self {
                res = T(
                    bs_non_central_chisquare_find_degreesOfFreedom_f(
                        Float(lambda),
                        Float(x),
                        Float(p))
                    )
            }
            else {
#if arch(x86_64) || arch(i386)
                res = T(
                    bs_non_central_chisquare_find_degreesOfFreedom_l(
                        Float80(lambda),
                        Float80(x),
                        Float80(p))
                    )
#else
                res = T(
                    bs_non_central_chisquare_find_degreesOfFreedom_d(
                        Double(lambda),
                        Double(x),
                        Double(p))
                    )
#endif
            }
            return Int(res)
        }

        public static func find_non_centrality(
            df: T,
            x: T,
            p: T) -> T {
            if T.self == Double.self {
                return T(
                    bs_non_central_chisquare_find_non_centrality_d(
                        Double(df),
                        Double(x),
                        Double(p))
                    )
            }
            else if T.self == Float.self {
                return T(
                    bs_non_central_chisquare_find_non_centrality_f(
                        Float(df),
                        Float(x),
                        Float(p))
                    )
            }
            else {
#if arch(x86_64) || arch(i386)
                return T(
                    bs_non_central_chisquare_find_non_centrality_l(
                        Float80(df),
                        Float80(x),
                        Float80(p))
                    )
#else
                return T(
                    bs_non_central_chisquare_find_non_centrality_d(
                        Double(df),
                        Double(x),
                        Double(p))
                    )
#endif
            }
        }

        
        // MARK: - Internal helpers

        @usableFromInline
        var moments: (mu1: T, mu2: T, mu3: T, mu4: T) {
            let v = self.degreesOfFreedom
            let lambda = self.nonCentrality

            let mu1 = v + lambda
            let mu2 = T(2) * (v + T(2) * lambda)
            let mu3 = T(8) * (v + T(3) * lambda)
            let base = v + T(2) * lambda
            let mu4 = T(48) * (v + T(4) * lambda) + T(12) * base * base
            return (mu1, mu2, mu3, mu4)
        }
    }
}
