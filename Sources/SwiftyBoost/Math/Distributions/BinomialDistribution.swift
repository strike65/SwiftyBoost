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

    /// A real-typed wrapper for the Binomial(n, p) distribution backed by Boost.Math.
    ///
    /// This type models the distribution of the number of successes in `nTrials` independent
    /// Bernoulli experiments where each trial has probability of success `pSuccess`.
    ///
    /// - Generic parameter T: The floating-point type used for evaluations and results.
    ///   Supported at runtime are:
    ///   - `Double` (preferred, highest precision on all Apple platforms),
    ///   - `Float`,
    ///   - `Float80` (on x86_64 only; falls back to `Double` bridge on other architectures).
    ///
    /// Conformance
    /// - `Sendable`: Instances are immutable and safe to use across concurrency domains.
    /// - `DistributionProtocol`: Exposes support bounds, density/mass (PDF), cumulative functions (CDF/SF),
    ///   quantiles, summary statistics, and hazard/CHF.
    ///
    /// Backing implementation
    /// - Core distribution behavior (PDF, CDF, quantiles, hazards, moments, etc.) is delegated to
    ///   an internal `Distribution.Dynamic<T>` which bridges to Boost.Math’s `binomial_distribution`.
    /// - Planning helpers (one-sided confidence bounds and trial-count solvers) call specific Boost
    ///   static utility functions via a C bridge (see bs_distribution_helpers.hxx).
    ///
    /// Domain and support
    /// - Support is the discrete lattice {0, 1, …, nTrials}. Evaluations at non-integer `x` follow
    ///   Boost’s continuous-extension semantics (internally using incomplete beta functions).
    ///
    /// Thread safety
    /// - This type is immutable and thread-safe. All methods are pure functions.
    public struct Binomial<T: Real & BinaryFloatingPoint & Sendable>: Sendable, DistributionProtocol {

        // MARK: - Types

        /// The real number type used by this binomial distribution.
        ///
        /// See the generic parameter `T` for supported concrete types.
        typealias RealType = T

        // MARK: - Parameters

        /// The total number of Bernoulli trials, `n ≥ 0`.
        ///
        /// - Note: Although `nTrials` is stored as `T`, it semantically represents a non-negative
        ///   integer count. Boost allows a continuous extension for some operations, but conceptually
        ///   this is an integer parameter.
        public let nTrials: T

        /// The probability of success for each trial, `0 ≤ p ≤ 1`.
        public let pSuccess: T

        /// Internal dynamic delegate that maps calls to Boost-backed implementations.
        private let dyn: Distribution.Dynamic<T>

        // MARK: - Initialization

        /// Creates a binomial distribution with the given number of trials and success probability.
        ///
        /// - Parameters:
        ///   - n: The number of Bernoulli trials (`n ≥ 0`).
        ///   - p: The probability of success for each trial (`0 ≤ p ≤ 1`).
        /// - Throws:
        ///   - `DistributionError.parameterNotPositive(name:value:)` if `n < 0`.
        ///   - `DistributionError.parameterOutOfRange(name:min:max:)` if `p` is outside `[0, 1]`.
        /// - Discussion:
        ///   The initializer validates parameters locally and then constructs a dynamic
        ///   Boost-backed delegate. All subsequent evaluations are forwarded to that delegate.
        public init(numberOfTrials n: T, probabibilityOfSuccess p: T) throws {
            guard n >= 0 else {
                throw DistributionError.parameterNotPositive(
                    name: "numberOfTrials", value: n
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

        // MARK: - DistributionProtocol — Support

        /// The lower bound of the support.
        ///
        /// For the binomial distribution this is typically 0.
        public var supportLowerBound: T { dyn.supportLowerBound }

        /// The upper bound of the support.
        ///
        /// For the binomial distribution this is typically `nTrials`.
        public var supportUpperBound: T { dyn.supportUpperBound }

        /// A convenient tuple of `(supportLowerBound, supportUpperBound)`.
        public var range: (lower: T, upper: T) { dyn.range }

        // MARK: - PDF / CDF / SF / Quantiles

        /// Probability mass (density) function evaluated at `x`.
        ///
        /// - Parameter x: Evaluation point. For integer `x` in `[0, nTrials]`,
        ///   this is the probability of exactly `x` successes.
        /// - Returns: The mass/density at `x`.
        /// - Throws: If `x` is outside the supported domain or if the underlying computation fails.
        public func pdf(_ x: T) throws -> T { try dyn.pdf(x) }

        /// Log of the probability mass (density) function at `x`.
        ///
        /// - Parameter x: Evaluation point.
        /// - Returns: `log(pdf(x))`, computed in a numerically stable fashion by the backend.
        /// - Throws: If `x` is outside the supported domain or if the underlying computation fails.
        public func logPdf(_ x: T) throws -> T { try dyn.logPdf(x) }

        /// Cumulative distribution function F(x) = P(X ≤ x).
        ///
        /// - Parameter x: Evaluation point.
        /// - Returns: A probability in [0, 1].
        /// - Throws: If `x` is outside the supported domain or if the underlying computation fails.
        public func cdf(_ x: T) throws -> T { try dyn.cdf(x) }

        /// Survival function S(x) = P(X > x) = 1 − F(x).
        ///
        /// - Parameter x: Evaluation point.
        /// - Returns: A probability in [0, 1].
        /// - Throws: If `x` is outside the supported domain or if the underlying computation fails.
        public func sf(_ x: T) throws -> T { try dyn.sf(x) }

        /// Lower-tail quantile function (inverse CDF).
        ///
        /// - Parameter p: A probability in [0, 1].
        /// - Returns: `x` such that `cdf(x) = p`.
        /// - Throws: If `p` is outside [0, 1] or if the computation fails to converge.
        public func quantile(_ p: T) throws -> T { try dyn.quantile(p) }

        /// Upper-tail quantile function (inverse survival).
        ///
        /// - Parameter q: An upper-tail probability in [0, 1].
        /// - Returns: `x` such that `sf(x) = q`.
        /// - Throws: If `q` is outside [0, 1] or if the computation fails to converge.
        public func quantileComplement(_ q: T) throws -> T {
            try dyn.quantileComplement(q)
        }

        // MARK: - Moments

        /// The expected value (mean) of the distribution.
        ///
        /// - Note:
        ///   For a binomial distribution, the theoretical mean is `nTrials * pSuccess`.
        ///   If a mean is not defined under the current parameterization or backend policy,
        ///   this may be `nil`.
        public var mean: T? { nTrials > 1 ? T(0) : nil }

        /// The variance of the distribution.
        ///
        /// - Note:
        ///   For a binomial distribution, the theoretical variance is `nTrials * pSuccess * (1 - pSuccess)`.
        public var variance: T? { dyn.variance }

        /// The mode (most probable number of successes), if uniquely defined.
        ///
        /// - Note:
        ///   For binomial, a common closed-form is `floor(p * (n + 1))`.
        public var mode: T? { dyn.mode }

        /// The median of the distribution.
        ///
        /// - Note:
        ///   Boost chooses a conventional median for discrete distributions; see Boost.Math docs.
        public var median: T { dyn.median }

        /// The skewness of the distribution, if defined.
        ///
        /// - Note:
        ///   For binomial, skewness is `(1 − 2p) / sqrt(n p (1 − p))` when defined.
        public var skewness: T? { dyn.skewness }

        /// The kurtosis (Pearson) of the distribution, if defined.
        ///
        /// - Note:
        ///   For binomial, kurtosis is `3 − 6/n + 1/(n p (1 − p))` when defined.
        public var kurtosis: T? { dyn.kurtosis }

        /// The excess kurtosis of the distribution, if defined.
        ///
        /// - Note:
        ///   For binomial, excess kurtosis is `(1 − 6 p (1 − p)) / (n p (1 − p))` when defined.
        public var kurtosisExcess: T? { dyn.kurtosisExcess }

        // MARK: - Hazards

        /// The (instantaneous) hazard function h(x) where defined.
        ///
        /// - Parameter x: Evaluation point.
        /// - Returns: The hazard at `x`.
        /// - Throws: If `x` is outside the supported domain or if the underlying computation fails.
        public func hazard(_ x: T) throws -> T { try dyn.hazard(x) }

        /// The cumulative hazard function H(x) where defined.
        ///
        /// - Parameter x: Evaluation point.
        /// - Returns: The cumulative hazard at `x)`.
        /// - Throws: If `x` is outside the supported domain or if the underlying computation fails.
        public func chf(_ x: T) throws -> T { try dyn.chf(x) }

        // MARK: - Lattice/discrete-only properties

        /// The lattice step for discrete support.
        ///
        /// - Note: For strictly discrete distributions like binomial, the natural step is 1.
        ///   This implementation returns `nil` and defers to the dynamic backend where applicable.
        public var latticeStep: T? { nil }

        /// The lattice origin for discrete support.
        ///
        /// - Note: For binomial, a natural origin is 0.
        ///   This implementation returns `nil` and defers to the dynamic backend where applicable.
        public var latticeOrigin: T? { nil }

        /// The Shannon entropy H[X] where defined (nats).
        public var entropy: T? { dyn.entropy }

        // MARK: - Planning helpers (one-sided bounds and trial-count solvers)
        //
        // These helpers expose Boost.Math’s static utilities to:
        // - Compute one-sided confidence bounds on the success probability p based on observed counts.
        // - Compute minimum/maximum trial counts to achieve certain risk levels for observing a given number of successes.
        //
        // Notes:
        // - `useJeffreys`: When true, uses Jeffreys prior interval; otherwise uses the Clopper–Pearson exact interval.
        // - All functions validate the probability argument is in [0, 1] where applicable.
        // - For generic `T`, the call is dispatched to the appropriate C bridge based on the concrete type of `T`.
        // - On non–x86_64 architectures, `Float80` is bridged via the `Double` implementation.

        /// Computes a one-sided lower confidence bound on the success fraction `p` given `n` trials and `k` observed successes.
        ///
        /// - Parameters:
        ///   - n: Number of trials (non-negative integer).
        ///   - k: Number of observed successes (0…n).
        ///   - p0: One-sided significance level (α) as a probability in [0, 1].
        ///   - useJeffreys: If `true`, use Jeffreys prior interval; else use Clopper–Pearson exact interval.
        /// - Returns: Lower bound on `p` in [0, 1], or `nan` if `p0` is outside [0, 1].
        /// - Discussion:
        ///   This is a convenience wrapper over Boost’s
        ///   `binomial_distribution<>::find_lower_bound_on_p(trials, successes, alpha, interval_type)`.
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

        /// Computes a one-sided upper confidence bound on the success fraction `p` given `n` trials and `k` observed successes.
        ///
        /// - Parameters:
        ///   - n: Number of trials (non-negative integer).
        ///   - k: Number of observed successes (0…n).
        ///   - p0: One-sided significance level (α) as a probability in [0, 1].
        ///   - useJeffreys: If `true`, use Jeffreys prior interval; else use Clopper–Pearson exact interval.
        /// - Returns: Upper bound on `p` in [0, 1], or `nan` if `p0` is outside [0, 1].
        /// - Discussion:
        ///   This is a convenience wrapper over Boost’s
        ///   `binomial_distribution<>::find_upper_bound_on_p(trials, successes, alpha, interval_type)`.
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

        /// Computes the minimum number of trials `n` required to be at most `alpha` risk of observing fewer than `s` successes,
        /// given a proposed success fraction `p0`.
        ///
        /// - Parameters:
        ///   - s: Target number of successes (events) you want to see.
        ///   - p0: Proposed success fraction `p` in [0, 1].
        ///   - alpha: One-sided risk level (probability) in [0, 1].
        /// - Returns: The smallest integer `n` meeting the criterion.
        /// - Throws: `DistributionError.parameterOutOfRange(name:min:max:)` if `p0` is outside [0, 1].
        /// - Discussion:
        ///   Wraps Boost’s `binomial_distribution<>::find_minimum_number_of_trials(k, p, alpha)`.
        public static func findMinimumNumberOfTrials(successes s: Int,
                                                     proposedSuccessFraction p0: T,
                                                     alpha: T) throws -> Int {
            guard p0 >= 0.0, p0 <= 1.0 else {
                throw DistributionError
                    .parameterOutOfRange(name: "p0", min: 0, max: 1)
            }
            let res : Int
            if T.self == Double.self {
                res = Int(ceil(bs_binomial_find_minimum_number_of_trials_d(Double(s), Double(p0), Double(alpha))))
            }
            else if T.self == Float.self {
                res = Int(ceil(bs_binomial_find_minimum_number_of_trials_f(Float(s), Float(p0), Float(alpha))))
            }
            else {
                #if arch(i386) || arch(x86_64)
                res = Int(ceil(bs_binomial_find_minimum_number_of_trials_f(Float80(s), Float80(p0), Float80(alpha))))
                #else
                res = Int(ceil(bs_binomial_find_minimum_number_of_trials_d(Double(s), Double(p0), Double(alpha))))
                #endif
            }
            return res
        }

        /// Computes the maximum number of trials `n` allowed to remain at most `alpha` risk of observing at least `s` successes,
        /// given a proposed success fraction `p0`.
        ///
        /// - Parameters:
        ///   - s: Target number of successes (events).
        ///   - p0: Proposed success fraction `p` in [0, 1].
        ///   - alpha: One-sided risk level (probability) in [0, 1].
        /// - Returns: The largest integer `n` meeting the criterion.
        /// - Throws: `DistributionError.parameterOutOfRange(name:min:max:)` if `p0` is outside [0, 1].
        /// - Discussion:
        ///   Wraps Boost’s `binomial_distribution<>::find_maximum_number_of_trials(k, p, alpha)`.
        public static func findMaximumNumberOfTrials(successes s: Int,
                                                     proposedSuccessFraction p0: T,
                                                     alpha: T) throws -> Int {
            guard p0 >= 0.0, p0 <= 1.0 else {
                throw DistributionError
                    .parameterOutOfRange(name: "p0", min: 0, max: 1)
            }
            if T.self == Double.self {
                return Int(floor(bs_binomial_find_maximum_number_of_trials_d(Double(s), Double(p0), Double(alpha))))
            }
            else if T.self == Float.self {
                return Int(floor(bs_binomial_find_maximum_number_of_trials_f(Float(s), Float(p0), Float(alpha))))
            }
            else {
                #if arch(i386) || arch(x86_64)
                return Int(floor(bs_binomial_find_maximum_number_of_trials_f(Float80(s), Float80(p0), Float80(alpha))))
                #else
                return Int(floor(bs_binomial_find_maximum_number_of_trials_d(Double(s), Double(p0), Double(alpha))))
                #endif
            }
        }
    }
}
