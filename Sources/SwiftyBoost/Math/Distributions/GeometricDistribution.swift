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
    /// A real-typed wrapper for the Geometric(p) distribution backed by Boost.Math.
    ///
    /// This type models the distribution of the number of failures observed before the first success
    /// in a sequence of independent Bernoulli trials with success probability `pSuccess` per trial.
    ///
    /// Conformance
    /// - `Sendable`: Instances are immutable and safe to use across concurrency domains.
    /// - `DistributionProtocol`: Exposes support bounds, density/mass (PDF), cumulative functions (CDF/SF),
    ///   quantiles, summary statistics, and hazard/CHF.
    ///
    /// Backing implementation
    /// - Core distribution behavior (PDF, CDF, quantiles, hazards, moments, etc.) is delegated to
    ///   an internal `Distribution.Dynamic<T>` which bridges to Boost.Math’s `geometric_distribution`.
    ///
    /// Domain and support
    /// - Interprets the geometric distribution in the “failures before first success” convention.
    /// - The support is the non-negative axis [0, +∞) with a discrete lattice at non-negative integers.
    ///   Boost allows a continuous extension for some operations; this wrapper forwards to Boost’s behavior.
    ///
    /// Thread safety
    /// - This type is immutable and thread-safe. All methods are pure functions.
    public struct Geometric<T: Real & BinaryFloatingPoint & Sendable>: Sendable, DistributionProtocol {

        // MARK: - Types

        /// The real number type used by this geometric distribution.
        ///
        /// Concrete types supported at runtime:
        /// - `Double` (preferred on Apple platforms),
        /// - `Float`,
        /// - `Float80` (x86_64 only; falls back to `Double` bridge on other architectures).
        public typealias RealType = T

        // MARK: - Parameters

        /// The per-trial probability of success, constrained to `0 ≤ p ≤ 1`.
        public let pSuccess: T

        /// Internal dynamic delegate that maps calls to Boost-backed implementations.
        private let dyn: Distribution.Dynamic<T>

        // MARK: - Initialization

        /// Creates a geometric distribution with the given probability of success per trial.
        ///
        /// - Parameter p: Probability of success per trial (`0 ≤ p ≤ 1`).
        /// - Throws: `DistributionError.parameterOutOfRange(name:min:max:)` if `p` is outside `[0, 1]`.
        /// - Discussion:
        ///   The initializer validates `p` locally and then constructs a dynamic Boost-backed delegate.
        ///   All subsequent evaluations are forwarded to that delegate.
        public init(probabibilityOfSuccess p: T) throws {
            guard p <= 1 && p >= 0 else {
                throw DistributionError
                    .parameterOutOfRange(name: "p", min: 0, max: 1)
            }
            self.pSuccess = p
            self.dyn = try Distribution.Dynamic<T>(
                distributionName: "geometric",
                parameters: [
                    "p": p
                ]
            )
        }

        // MARK: - DistributionProtocol — Support

        /// The lower bound of the support (typically 0).
        public var supportLowerBound: T { dyn.supportLowerBound }

        /// The upper bound of the support (typically `+infinity` where supported).
        public var supportUpperBound: T { dyn.supportUpperBound }

        /// A convenient tuple of `(supportLowerBound, supportUpperBound)`.
        public var range: (lower: T, upper: T) { dyn.range }

        // MARK: - PDF / CDF / SF / Quantiles

        /// Probability mass (density) evaluated at `x`.
        ///
        /// - Parameter x: Evaluation point. For integer `x ≥ 0`, this is the probability of observing
        ///   exactly `x` failures before the first success.
        /// - Returns: The mass/density at `x`.
        /// - Throws: If `x` is outside the supported domain or if the underlying computation fails.
        public func pdf(_ x: T) throws -> T { try dyn.pdf(x) }

        /// Log of the probability mass (density) at `x`.
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

        /// The expected value (mean) of the distribution, if defined.
        ///
        /// - Note:
        ///   For a geometric distribution (failures before first success), the theoretical mean is `(1 − p) / p`.
        ///   This wrapper forwards to the Boost backend and returns `nil` if the backend reports it as undefined.
        public var mean: T? { dyn.mean }

        /// The variance of the distribution, if defined.
        ///
        /// - Note:
        ///   For geometric, the theoretical variance is `(1 − p) / p^2`.
        public var variance: T? { dyn.variance }

        /// The mode (most probable number of failures), if uniquely defined.
        ///
        /// - Note:
        ///   For geometric, the mode is typically 0 under this parameterization.
        public var mode: T? { dyn.mode }

        /// The median of the distribution.
        ///
        /// - Note:
        ///   Boost chooses a conventional median for discrete distributions; see Boost.Math docs.
        public var median: T { dyn.median }

        /// The skewness of the distribution, if defined.
        ///
        /// - Note:
        ///   For geometric, skewness is `(2 − p) / sqrt(1 − p)`.
        public var skewness: T? { dyn.skewness }

        /// The kurtosis (Pearson) of the distribution, if defined.
        ///
        /// - Note:
        ///   For geometric, kurtosis is `9 − 6p + p^2 / (1 − p)`.
        public var kurtosis: T? { dyn.kurtosis }

        /// The excess kurtosis of the distribution, if defined.
        ///
        /// - Note:
        ///   For geometric, excess kurtosis is `6 + p^2 / (1 − p)`.
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
        /// - Returns: The cumulative hazard at `x`.
        /// - Throws: If `x` is outside the supported domain or if the underlying computation fails.
        public func chf(_ x: T) throws -> T { try dyn.chf(x) }

        // MARK: - Lattice/discrete-only properties

        /// The lattice step for discrete support (spacing between adjacent support points).
        ///
        /// - Note: For strictly discrete distributions like geometric, the natural step is 1.
        ///   This implementation returns `nil` and defers to the dynamic backend where applicable.
        public var latticeStep: T? { nil }

        /// The lattice origin for discrete support.
        ///
        /// - Note: For geometric, a natural origin is 0.
        ///   This implementation returns `nil` and defers to the dynamic backend where applicable.
        public var latticeOrigin: T? { nil }

        /// The Shannon entropy H[X] of the geometric distribution (in nats).
        ///
        /// - Definition:
        ///   For the “failures before first success” parameterization with success probability `p = pSuccess`,
        ///   the entropy is
        ///   H = −ln(p) − ((1 − p) / p) · ln(1 − p),
        ///   where all logarithms are natural logs (base e).
        ///
        /// - Units:
        ///   Returned in nats. To convert to bits, divide by ln(2).
        ///
        /// - Edge cases:
        ///   - If `p == 1`, the distribution is degenerate at 0 failures and the entropy is 0.
        ///   - If `0 < p < 1`, this returns a finite positive value.
        ///   - If `p == 0`, the entropy tends to +infinity; this will surface as `.infinity`.
        public var entropy: T? {
            if self.pSuccess == 1 { return 0 }
            if self.pSuccess.isZero { return .infinity }
            let ln1mp = T.log(onePlus: -self.pSuccess)
            let hNats = -T.log(self.pSuccess) - ((1 - self.pSuccess ) / self.pSuccess) * ln1mp
            return hNats
        }
        /// Indicates whether this distribution is discrete (`true`) or continuous (`false`).
        ///
        /// Geometric distributions are discrete on the non-negative integers, so this always returns `true`.
        public var isDiscrete: Bool { dyn.isDiscrete }

        /// Computes the Kullback–Leibler divergence `D_KL(self || other)` when defined.
        ///
        /// - Parameters:
        ///   - other: The reference geometric distribution *Q*.
        ///   - options: Summation/integration configuration; defaults to ``Distribution/KLDivergenceOptions/automatic()``.
        /// - Returns: The divergence in nats, or `nil` if it cannot be evaluated.
        /// - Throws: Rethrows any backend or quadrature errors.
        public func klDivergence<D>(
            relativeTo other: D,
            options: Distribution.KLDivergenceOptions<T>
        ) throws -> T? where D: DistributionProtocol, D.RealType == T {
            try DistributionKLDivergenceHelper.evaluate(lhs: self, rhs: other, options: options)
        }
// MARK: - Planning helpers (one-sided bounds and trial-count solvers)
        //
        // These helpers expose Boost.Math’s static utilities to:
        // - Compute one-sided confidence bounds on the success probability p based on observed counts.
        // - Compute minimum/maximum trial counts to achieve certain risk levels for observing a given number of failures.
        //
        // Notes:
        // - All functions validate the probability argument is in [0, 1] where applicable.
        // - For generic `T`, the call is dispatched to the appropriate C bridge based on the concrete type of `T`.
        // - On non–x86_64 architectures, `Float80` is bridged via the `Double` implementation.
        // - These helpers return `nan` if the probability argument is outside [0, 1] rather than throwing.

        /// Computes a one-sided lower confidence bound on the success fraction `p` given `n` observed trials.
        ///
        /// - Parameters:
        ///   - n: Number of trials (non-negative integer), interpreted per Boost’s geometric utility.
        ///   - alpha: One-sided significance level (α) as a probability in [0, 1].
        /// - Returns: Lower bound on `p` in [0, 1], or `nan` if `alpha` is outside [0, 1].
        /// - Discussion:
        ///   Wraps Boost’s `geometric_distribution<>::find_lower_bound_on_p(trials, alpha)`.
        public static func findLowerBoundOnP(nTrials n: Int,
                                             alpha: T) -> T {
            guard alpha >= 0.0, alpha <= 1.0 else { return .nan  }
            if T.self == Double.self {
               return T(bs_geometric_find_lower_bound_on_p_d(
                   Double(n),
                   Double(alpha)))
            }
            else if T.self == Float.self {
                return T(bs_geometric_find_lower_bound_on_p_f(
                    Float(n),
                    Float(alpha)))
            }
            else {
                #if arch(i386) || arch(x86_64)
                return T(bs_geometric_find_lower_bound_on_p_l(
                    Float80(n),
                    Float80(alpha)))
                #else
                return T(bs_geometric_find_lower_bound_on_p_d(
                    Double(n),
                    Double(alpha)))
                #endif
            }
        }

        /// Computes a one-sided upper confidence bound on the success fraction `p` given `n` observed trials.
        ///
        /// - Parameters:
        ///   - n: Number of trials (non-negative integer), interpreted per Boost’s geometric utility.
        ///   - alpha: One-sided significance level (α) as a probability in [0, 1].
        /// - Returns: Upper bound on `p` in [0, 1], or `nan` if `alpha` is outside [0, 1].
        /// - Discussion:
        ///   Wraps Boost’s `geometric_distribution<>::find_upper_bound_on_p(trials, alpha)`.
        public static func findUpperBoundOnP(nTrials n: Int,
                                             alpha: T) -> T {
            guard alpha >= 0.0, alpha <= 1.0 else { return .nan }
            if T.self == Double.self {
               return T(bs_geometric_find_upper_bound_on_p_d(
                   Double(n),
                   Double(alpha)))
            }
            else if T.self == Float.self {
                return T(bs_geometric_find_upper_bound_on_p_f(
                    Float(n),
                    Float(alpha)))
            }
            else {
                #if arch(i386) || arch(x86_64)
                return T(bs_geometric_find_upper_bound_on_p_l(
                    Float80(n),
                    Float80(alpha)))
                #else
                return T(bs_geometric_find_upper_bound_on_p_d(
                    Double(n),
                    Double(alpha)))
                #endif
            }
        }

        /// Computes the minimum number of trials `n` required to be at most `alpha` risk of observing more than `f` failures
        /// before the first success, given a proposed success fraction `p0`.
        ///
        /// - Parameters:
        ///   - f: Target number of failures you want to stay at or below (before first success).
        ///   - p0: Proposed success fraction `p` in [0, 1].
        ///   - alpha: One-sided risk level (probability) in [0, 1].
        /// - Returns: The smallest integer `n` meeting the criterion.
        /// - Throws:
        ///   - `DistributionError.parameterOutOfRange(name:min:max:)` if `p0` is outside [0, 1].
        ///   - `DistributionError.parameterOutOfRange(name:min:max:)` if `alpha` is outside [0, 1].
        /// - Discussion:
        ///   Wraps Boost’s `geometric_distribution<>::find_minimum_number_of_trials(failures, p, alpha)`.
        public static func findMinimumNumberOfTrials(failures f: Int,
                                                     proposedSuccessFraction p0: T,
                                                     alpha: T) throws -> Int {
            guard p0 >= 0.0, p0 <= 1.0 else {
                throw DistributionError
                    .parameterOutOfRange(name: "p0", min: 0, max: 1)
            }
            guard alpha >= 0.0, alpha <= 1.0 else {
                throw DistributionError
                    .parameterOutOfRange(name: "alpha", min: 0, max: 1)
            }
            let res : Int
            if T.self == Double.self {
                res = Int(ceil(bs_geometric_find_minimum_number_of_trials_d(Double(f), Double(p0), Double(alpha))))
            }
            else if T.self == Float.self {
                res = Int(ceil(bs_geometric_find_minimum_number_of_trials_f(Float(f), Float(p0), Float(alpha))))
            }
            else {
                #if arch(i386) || arch(x86_64)
                res = Int(ceil(bs_geometric_find_minimum_number_of_trials_l(Float80(f), Float80(p0), Float80(alpha))))
                #else
                res = Int(ceil(bs_geometric_find_minimum_number_of_trials_d(Double(f), Double(p0), Double(alpha))))
                #endif
            }
            return res
        }

        /// Computes the maximum number of trials `n` allowed to remain at most `alpha` risk of observing at least `f` failures
        /// before the first success, given a proposed success fraction `p0`.
        ///
        /// - Parameters:
        ///   - f: Target number of failures (before first success).
        ///   - p0: Proposed success fraction `p` in [0, 1].
        ///   - alpha: One-sided risk level (probability) in [0, 1].
        /// - Returns: The largest integer `n` meeting the criterion.
        /// - Throws:
        ///   - `DistributionError.parameterOutOfRange(name:min:max:)` if `p0` is outside [0, 1].
        ///   - `DistributionError.parameterOutOfRange(name:min:max:)` if `alpha` is outside [0, 1].
        /// - Discussion:
        ///   Wraps Boost’s `geometric_distribution<>::find_maximum_number_of_trials(failures, p, alpha)`.
        public static func findMaximumNumberOfTrials(failures f: Int,
                                                     proposedSuccessFraction p0: T,
                                                     alpha: T) throws -> Int {
            guard p0 >= 0.0, p0 <= 1.0 else {
                throw DistributionError
                    .parameterOutOfRange(name: "p0", min: 0, max: 1)
            }
            guard alpha >= 0.0, alpha <= 1.0 else {
                throw DistributionError
                    .parameterOutOfRange(name: "alpha", min: 0, max: 1)
            }
            let res : Int
            if T.self == Double.self {
                res = Int(floor(bs_geometric_find_maximum_number_of_trials_d(Double(f), Double(p0), Double(alpha))))
            }
            else if T.self == Float.self {
                res = Int(floor(bs_geometric_find_maximum_number_of_trials_f(Float(f), Float(p0), Float(alpha))))
            }
            else {
                #if arch(i386) || arch(x86_64)
                res = Int(floor(bs_geometric_find_maximum_number_of_trials_l(Float80(f), Float80(p0), Float80(alpha))))
                #else
                res = Int(floor(bs_geometric_find_maximum_number_of_trials_d(Double(f), Double(p0), Double(alpha))))
                #endif
            }
            return res
        }
    }
}
