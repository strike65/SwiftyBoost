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
    public struct Binomial<T: Real & BinaryFloatingPoint & Sendable>: Sendable,
        DistributionProtocol
    {
        typealias RealType = T
        public let nTrials: T
        public let pSuccess: T
        private let dyn: Distribution.Dynamic<T>

        public init(numberOfTrials n: T, probabibilityOfSuccess p: T) throws {
            guard n >= 0 else {
                throw DistributionError.parameterNotPositive(
                    name: "degreesOfFreedom"
                )
            }
            guard p <= 1 && p >= 0 else {
                throw DistributionError
                    .parameterOutOfRange(name: "p", min: 0, max: 1)
            }
            self.pSuccess = p
            self.nTrials = n
            self.dyn = try Distribution.Dynamic<T>(
                distributionName: "binomial",
                parameters: [
                    "n": n,
                    "p": p
                ]
            )
        }

        // MARK: DistributionProtocol — Support
        public var supportLowerBound: T { dyn.supportLowerBound }
        public var supportUpperBound: T { dyn.supportUpperBound }
        public var range: (lower: T, upper: T) { dyn.range }

        // MARK: PDF/CDF/SF/Quantile
        public func pdf(_ x: T) throws -> T { try dyn.pdf(x) }
        public func logPdf(_ x: T) throws -> T { try dyn.logPdf(x) }
        public func cdf(_ x: T) throws -> T { try dyn.cdf(x) }
        public func sf(_ x: T) throws -> T { try dyn.sf(x) }
        public func quantile(_ p: T) throws -> T { try dyn.quantile(p) }
        public func quantileComplement(_ q: T) throws -> T {
            try dyn.quantileComplement(q)
        }

        // MARK: Moments
        public var mean: T? { nTrials > 1 ? T(0) : nil }
        public var variance: T? { dyn.variance }
        public var mode: T? { dyn.mode }
        public var median: T { dyn.median }
        public var skewness: T? { dyn.skewness }
        public var kurtosis: T? { dyn.kurtosis }
        public var kurtosisExcess: T? { dyn.kurtosisExcess }

        // MARK: Hazards
        public func hazard(_ x: T) throws -> T { try dyn.hazard(x) }
        public func chf(_ x: T) throws -> T { try dyn.chf(x) }

        // Lattice/discrete-only properties (continuous ⇒ nil)
        public var latticeStep: T? { nil }
        public var latticeOrigin: T? { nil }
        public var entropy: T? { dyn.entropy }

        // MARK: Planning helper
        /// Find degrees of freedom ν given effect size, α, β, and σ.
        public static func findLowerBoundOnP(nTrials n: Int,
                                             nSuccesses k: Int,
                                             proposedSuccessFraction p0: T,
                                             useJeffreys: Bool = false) -> T {
            guard p0 >= 0.0, p0 <= 1.0 else { return .nan }
            if T.self == Double.self {
               return T(bs_binomial_find_lower_bound_on_p_d(
                   Double(n),
                   Double(k),
                   Double(p0),
                   useJeffreys))
            }
            else if T.self == Float.self {
                return T(bs_binomial_find_lower_bound_on_p_f(
                    Float(n),
                    Float(k),
                    Float(p0),
                    useJeffreys))
            }
            else {
                #if arch(i386) || arch(x86_64)
                return T(bs_binomial_find_lower_bound_on_p_f(
                    Float80(n),
                    Float80(k),
                    Float80(p0),
                    useJeffreys))
                #else
                return T(bs_binomial_find_lower_bound_on_p_d(
                    Double(n),
                    Double(k),
                    Double(p0),
                    useJeffreys))
                #endif
            }
        }

        public static func findUpperBoundOnP(nTrials n: Int,
                                             nSuccesses k: Int,
                                             proposedSuccessFraction p0: T,
                                             useJeffreys: Bool = false) -> T {
            guard p0 >= 0.0, p0 <= 1.0 else { return .nan }
            if T.self == Double.self {
               return T(bs_binomial_find_upper_bound_on_p_d(
                   Double(n),
                   Double(k),
                   Double(p0),
                   useJeffreys))
            }
            else if T.self == Float.self {
                return T(bs_binomial_find_upper_bound_on_p_f(
                    Float(n),
                    Float(k),
                    Float(p0),
                    useJeffreys))
            }
            else {
                #if arch(i386) || arch(x86_64)
                return T(bs_binomial_find_upper_bound_on_p_f(
                    Float80(n),
                    Float80(k),
                    Float80(p0),
                    useJeffreys))
                #else
                return T(bs_binomial_find_upper_bound_on_p_d(
                    Double(n),
                    Double(k),
                    Double(p0),
                    useJeffreys))
                #endif
            }
        }

        public static func findMinimumNumberOfTrials(successes s: Int,
                                                     proposedSuccessFraction p0: T,
                                                     alpha: T) -> T {
            guard p0 >= 0.0, p0 <= 1.0 else { return .nan }
            if T.self == Double.self {
                return T(bs_binomial_find_minimum_number_of_trials_d(Double(s), Double(p0), Double(alpha)))
            }
            else if T.self == Float.self {
                return T(bs_binomial_find_minimum_number_of_trials_f(Float(s), Float(p0), Float(alpha)))
            }
            else {
                #if arch(i386) || arch(x86_64)
                return T(bs_binomial_find_minimum_number_of_trials_f(Float80(s), Float80(p0), Float80(alpha)))
                #else
                return T(bs_binomial_find_minimum_number_of_trials_d(Double(s), Double(p0), Double(alpha)))
                #endif
            }
        }

        public static func findMaximumNumberOfTrials(successes s: Int,
                                                     proposedSuccessFraction p0: T,
                                                     alpha: T) -> T {
            guard p0 >= 0.0, p0 <= 1.0 else { return .nan }
            if T.self == Double.self {
                return T(bs_binomial_find_maximum_number_of_trials_d(Double(s), Double(p0), Double(alpha)))
            }
            else if T.self == Float.self {
                return T(bs_binomial_find_maximum_number_of_trials_f(Float(s), Float(p0), Float(alpha)))
            }
            else {
                #if arch(i386) || arch(x86_64)
                return T(bs_binomial_find_maximum_number_of_trials_f(Float80(s), Float80(p0), Float80(alpha)))
                #else
                return T(bs_binomial_find_maximum_number_of_trials_d(Double(s), Double(p0), Double(alpha)))
                #endif
            }
        }
    }
}
