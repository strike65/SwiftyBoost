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
    /// Negative binomial distribution (failures before a fixed number of successes).
    ///
    /// The parameters follow the standard “Pascal” convention:
    /// - `successes`: the (possibly non-integer) target number of successes `r > 0`.
    /// - `probabilityOfSuccess`: the probability of success per trial `p`, with `0 < p ≤ 1`.
    ///
    /// The support is the non-negative integers counting the number of failures observed before
    /// the `successes`-th success occurs.
    public struct NegativeBinomial<T: Real & BinaryFloatingPoint & Sendable>: Sendable, DistributionProtocol {
        typealias RealType = T

        /// Target number of successes `r`.
        public let successes: T

        /// Probability of success per trial `p`.
        public let probabilityOfSuccess: T

        private let dyn: Distribution.Dynamic<T>

        /// Creates a negative binomial distribution parameterised by successes and success probability.
        ///
        /// - Parameters:
        ///   - successes: Target number of successes `r`; must be strictly positive.
        ///   - probabilityOfSuccess: Success probability per trial `p`; must satisfy `0 < p ≤ 1`.
        /// - Throws: `DistributionError.parameterNotPositive` if `successes <= 0`,
        ///           `DistributionError.parameterOutOfRange` if `probabilityOfSuccess` is outside `(0, 1]`.
        public init(successes r: T, probabilityOfSuccess p: T) throws {
            guard r > 0 else {
                throw DistributionError.parameterNotPositive(name: "successes", value: r)
            }
            guard p > 0 && p <= 1 else {
                throw DistributionError.parameterOutOfRange(name: "probabilityOfSuccess", min: 0, max: 1)
            }
            self.successes = r
            self.probabilityOfSuccess = p
            self.dyn = try Distribution.Dynamic<T>(
                distributionName: "negative_binomial",
                parameters: [
                    "successes": r,
                    "p": p
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

        // MARK: Moments and summary statistics

        public var mean: T? { dyn.mean }
        public var variance: T? { dyn.variance }
        public var mode: T? { dyn.mode }
        public var median: T { dyn.median }
        public var skewness: T? { dyn.skewness }
        public var kurtosis: T? { dyn.kurtosis }
        public var kurtosisExcess: T? { dyn.kurtosisExcess }
        public var entropy: T? {
            let rD = self.successes
            let pD = self.probabilityOfSuccess

            guard pD > 0, pD <= 1, rD > 0,
                  pD.isFinite, rD.isFinite else {
                return nil
            }

            if pD == 1 {
                return T.zero
            }

            let qD = 1 - pD

            var logProbability = rD * T.log(pD)  // log P(K = 0)
            var probability: T = T.exp(logProbability)
            if !probability.isFinite {
                return nil
            }

            var entropyAcc:T = 0.0
            var correction:T = 0.0
            var cumulative:T = probability
            var iterations = 0

            let tailTolerance:T = 1e-12
            let termTolerance: T = 1e-20
            let maxIterations = 1_000_000
            let pTermTol: T = 1e-9
            while true {
                if probability > 0 {
                    let term = -probability * T.log(probability)
                    let sum = entropyAcc + term
                    if abs(entropyAcc) >= abs(term) {
                        correction += (entropyAcc - sum) + term
                    } else {
                        correction += (term - sum) + entropyAcc
                    }
                    entropyAcc = sum
                }

                let tail = max(0.0, 1.0 - cumulative)
                if tail <= tailTolerance || !tail.isFinite {
                    break
                }

                iterations += 1
                if iterations >= maxIterations {
                    break
                }

                logProbability += T.log(qD)
                logProbability += T.log(rD + T(iterations - 1))
                logProbability -= T.log(T(iterations))
                probability = T.exp(logProbability)

                if probability <= 0 || !probability.isFinite {
                    break
                }

                let nextCumulative = cumulative + probability
                if nextCumulative == cumulative {
                    break
                }
                cumulative = min(1.0, nextCumulative)

                if probability <= termTolerance && tail <= pTermTol {
                    break
                }
            }

            let result = entropyAcc + correction
            if result.isFinite {
                return result
            } else {
                return nil
            }
        }

        // MARK: Lattice metadata

        public var latticeStep: T? { dyn.latticeStep }
        public var latticeOrigin: T? { dyn.latticeOrigin }
        public var isDiscrete: Bool { dyn.isDiscrete }

        /// Delegates KL divergence to the dynamic backend (numerical integration).
        public func klDivergence(
            relativeTo other: Self,
            options: Distribution.KLDivergenceOptions<T>
        ) throws -> T? {
            try dyn.klDivergence(relativeTo: other.dyn, options: options)
        }
        

        public static func findLowerBoundOnP(nTrials n: Int,
                                             nSuccesses k: Int,
                                             proposedSuccessFraction p0: T) -> T {
            guard p0 >= 0.0, p0 <= 1.0 else { return .nan }
            if T.self == Double.self {
               return T(bs_negative_binomial_find_lower_bound_on_p_d(
                   Double(n),
                   Double(k),
                   Double(p0)))
            }
            else if T.self == Float.self {
                return T(bs_negative_binomial_find_lower_bound_on_p_f(
                    Float(n),
                    Float(k),
                    Float(p0)))
            }
            else {
                #if arch(i386) || arch(x86_64)
                return T(bs_negative_binomial_find_lower_bound_on_p_l(
                    Float80(n),
                    Float80(k),
                    Float80(p0)))
                #else
                return T(bs_negative_binomial_find_lower_bound_on_p_d(
                    Double(n),
                    Double(k),
                    Double(p0)))
                #endif
            }
        }
        public static func findUpperBoundOnP(nTrials n: Int,
                                             nSuccesses k: Int,
                                             proposedSuccessFraction p0: T) -> T {
            guard p0 >= 0.0, p0 <= 1.0 else { return .nan }
            if T.self == Double.self {
               return T(bs_negative_binomial_find_upper_bound_on_p_d(
                   Double(n),
                   Double(k),
                   Double(p0)))
            }
            else if T.self == Float.self {
                return T(bs_negative_binomial_find_upper_bound_on_p_f(
                    Float(n),
                    Float(k),
                    Float(p0)))
            }
            else {
                #if arch(i386) || arch(x86_64)
                return T(bs_negative_binomial_find_upper_bound_on_p_l(
                    Float80(n),
                    Float80(k),
                    Float80(p0)))
                #else
                return T(bs_negative_binomial_find_upper_bound_on_p_d(
                    Double(n),
                    Double(k),
                    Double(p0)))
                #endif
            }
        }
        
        public static func findMinimumNumberOfTrials(successes s: Int,
                                                     proposedSuccessFraction p0: T,
                                                     alpha: T) throws -> Int {
            guard p0 >= 0.0, p0 <= 1.0 else {
                throw DistributionError
                    .parameterOutOfRange(name: "p0", min: 0, max: 1)
            }
            let res : Int
            if T.self == Double.self {
                res = Int(ceil(bs_negative_binomial_find_minimum_number_of_trials_d(Double(s), Double(p0), Double(alpha))))
            }
            else if T.self == Float.self {
                res = Int(ceil(bs_negative_binomial_find_minimum_number_of_trials_f(Float(s), Float(p0), Float(alpha))))
            }
            else {
                #if arch(i386) || arch(x86_64)
                res = Int(ceil(bs_negative_binomial_find_minimum_number_of_trials_l(Float80(s), Float80(p0), Float80(alpha))))
                #else
                res = Int(ceil(bs_negative_binomial_find_minimum_number_of_trials_d(Double(s), Double(p0), Double(alpha))))
                #endif
            }
            return res
        }

        public static func findMaximumNumberOfTrials(successes s: Int,
                                                     proposedSuccessFraction p0: T,
                                                     alpha: T) throws -> Int {
            guard p0 >= 0.0, p0 <= 1.0 else {
                throw DistributionError
                    .parameterOutOfRange(name: "p0", min: 0, max: 1)
            }
            let res : Int
            if T.self == Double.self {
                res = Int(floor(bs_negative_binomial_find_maximum_number_of_trials_d(Double(s), Double(p0), Double(alpha))))
            }
            else if T.self == Float.self {
                res = Int(floor(bs_negative_binomial_find_maximum_number_of_trials_f(Float(s), Float(p0), Float(alpha))))
            }
            else {
                #if arch(i386) || arch(x86_64)
                res = Int(floor(bs_negative_binomial_find_maximum_number_of_trials_l(Float80(s), Float80(p0), Float80(alpha))))
                #else
                res = Int(floor(bs_negative_binomial_find_maximum_number_of_trials_d(Double(s), Double(p0), Double(alpha))))
                #endif
            }
            return res
        }


    }
}
