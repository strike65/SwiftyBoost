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
//  BetaDistribution.swift
//
//  Overview
//  --------
//  A strongly-typed Swift wrapper around Boost.Math’s beta distribution, bridged
//  through CBoostBridge. This type exposes a uniform API (PDF, CDF, quantiles,
//  summary statistics, hazard functions, and entropy) via a generic interface
//  that works with Float, Double, and (on x86) Float80.
//
//  Usage
//  -----
//  do {
//      let dist = try Distribution.Beta<Double>(alpha: 2, beta: 5)
//      let p = try dist.cdf(0.3)
//      let x = try dist.quantile(0.95)
//      let m = dist.mean            // ≈ 0.2857
//      let v = dist.variance        // ≈ 0.02551
//  } catch {
//      // Handle DistributionError
//  }
//
//  Notes
//  -----
//  - Parameter validation follows Boost’s constraints: alpha > 0, beta > 0, both finite.
//  - The support is [0, 1]. Values outside this range yield zero density and CDF limits.
//  - Computations are delegated to Distribution.Dynamic<T>, which binds to Boost.
//  - Instances are immutable and Sendable; they can be shared across tasks safely.
//
//  Thread-safety
//  -------------
//  Instances are immutable. Calls into the underlying dynamic distribution wrapper are
//  expected to be thread-safe when used concurrently across independent instances.
//
import CBoostBridge
import Foundation

public extension Distribution {
    /// The Beta distribution on the unit interval [0, 1].
    ///
    /// Definition
    /// - Support: x ∈ [0, 1]
    /// - Parameters: α > 0 (alpha), β > 0 (beta)
    /// - PDF:
    ///   f(x; α, β) = x^(α − 1) (1 − x)^(β − 1) / B(α, β), for x ∈ (0, 1),
    ///   where B(α, β) is the beta function.
    ///
    /// Common uses
    /// - Bayesian modeling of probabilities and proportions.
    /// - Conjugate prior for Bernoulli, binomial, and related likelihoods.
    ///
    /// Generic type
    /// - T: A BinaryFloatingPoint & Sendable type (Float, Double, and on x86, Float80).
    ///
    /// Bridging
    /// - All operations delegate to a dynamically constructed Boost beta distribution via
    ///   Distribution.Dynamic<T> (CBoostBridge). This ensures consistent numerical behavior
    ///   and error handling with Boost.Math.
    ///
    /// Error handling
    /// - Throws DistributionError for invalid parameters or domain violations in function calls.
    /// - Summary statistics that are undefined or non-finite for the chosen T may return nil.
    ///
    /// Concurrency
    /// - Immutable and Sendable. Safe to share across tasks/threads.
    struct Beta<T: BinaryFloatingPoint & Sendable>: Sendable, DistributionProtocol {
        /// The floating-point type used by this distribution.
        typealias Real = T

        /// The first (α) shape parameter. Must be strictly greater than zero and finite.
        public let alpha: T
        /// The second (β) shape parameter. Must be strictly greater than zero and finite.
        public let beta: T
        /// Internal dynamic distribution wrapper bound to Boost's beta distribution.
        private let dyn: Distribution.Dynamic<T>

        // MARK: - Initialization

        /// Creates a Beta(α, β) distribution with the given shape parameters.
        ///
        /// - Parameters:
        ///   - alpha: Shape parameter α. Must be strictly positive and finite.
        ///   - beta: Shape parameter β. Must be strictly positive and finite.
        ///
        /// - Throws:
        ///   - DistributionError.invalidCombination if alpha <= 0 or beta <= 0.
        ///   - DistributionError.parameterNotFinite if alpha or beta is not finite.
        ///   - DistributionError.generalError if the underlying distribution could not be initialized.
        ///
        /// Notes:
        /// - The support is the closed interval [0, 1].
        /// - For extreme values of α or β, numerical stability is delegated to Boost.
        public init(alpha: T = 0, beta: T = 0) throws {
            guard alpha > 0 else { throw DistributionError.invalidCombination(message: "alpha must be > 0") }
            guard beta > 0 else { throw DistributionError.invalidCombination(message: "beta must be > 0") }
            guard alpha.isFinite else { throw DistributionError.parameterNotFinite(name: "alpha") }
            guard beta.isFinite else { throw DistributionError.parameterNotFinite(name: "beta") }
            self.alpha = alpha
            self.beta = beta
            self.dyn = try Distribution.Dynamic<T>(
                distributionName: "beta",
                parameters: ["alpha": alpha, "beta": beta]
            )
        }

        // MARK: - Support

        /// The lower bound of the support.
        ///
        /// For the Beta distribution this is 0.
        public var supportLowerBound: T { dyn.supportLowerBound }

        /// The upper bound of the support.
        ///
        /// For the Beta distribution this is 1.
        public var supportUpperBound: T { dyn.supportUpperBound }

        /// The support as a `(lower, upper)` tuple.
        ///
        /// For the Beta distribution this should be `(0, 1)`.
        public var range: (lower: T, upper: T) { dyn.range }

        // MARK: - Core functions

        /// The probability density function (PDF) at `x`.
        ///
        /// - Parameter x: Point of evaluation. Values outside `[0, 1]` should return zero.
        /// - Returns: f(x; α, β).
        /// - Throws: DistributionError from the underlying implementation if evaluation fails.
        public func pdf(_ x: T) throws -> T { try dyn.pdf(x) }

        /// The natural logarithm of the PDF at `x`.
        ///
        /// - Parameter x: Point of evaluation.
        /// - Returns: log f(x; α, β).
        /// - Throws: DistributionError from the underlying implementation if evaluation fails.
        public func logPdf(_ x: T) throws -> T { try dyn.logPdf(x) }

        /// The cumulative distribution function (CDF).
        ///
        /// - Parameter x: Point of evaluation.
        /// - Returns: P(X ≤ x).
        /// - Throws: DistributionError from the underlying implementation if evaluation fails.
        public func cdf(_ x: T) throws -> T { try dyn.cdf(x) }

        /// The survival function (SF), i.e. the complementary CDF.
        ///
        /// - Parameter x: Point of evaluation.
        /// - Returns: P(X > x) = 1 − F(x).
        /// - Throws: DistributionError from the underlying implementation if evaluation fails.
        public func sf(_ x: T) throws -> T { try dyn.sf(x) }

        /// The lower-tail quantile function (inverse CDF).
        ///
        /// - Parameter p: A probability in [0, 1].
        /// - Returns: x such that P(X ≤ x) = p.
        /// - Throws: DistributionError from the underlying implementation if evaluation fails.
        public func quantile(_ p: T) throws -> T { try dyn.quantile(p) }

        /// The upper-tail quantile function (inverse survival function).
        ///
        /// - Parameter q: A probability in [0, 1].
        /// - Returns: x such that P(X > x) = q.
        /// - Throws: DistributionError from the underlying implementation if evaluation fails.
        public func quantileComplement(_ q: T) throws -> T { try dyn.quantileComplement(q) }

        // MARK: - Summary statistics

        /// The mean of the distribution.
        ///
        /// For α, β > 0, mean = α / (α + β).
        ///
        /// - Note: Returns nil if not finite for the chosen numeric type.
        public var mean: T? { dyn.mean }

        /// The variance of the distribution.
        ///
        /// For α, β > 0, variance = αβ / ((α + β)^2 (α + β + 1)).
        ///
        /// - Note: Returns nil if not finite for the chosen numeric type.
        public var variance: T? { dyn.variance }

        /// The mode of the distribution.
        ///
        /// For α > 1 and β > 1, mode = (α − 1) / (α + β − 2). Otherwise, modes occur at the boundary.
        ///
        /// - Note: The exact behavior reflects the underlying Boost implementation.
        public var mode: T? { dyn.mode }

        /// The median of the distribution.
        ///
        /// There is no general closed form; obtained via the underlying implementation.
        public var median: T { dyn.median }

        /// The skewness of the distribution.
        ///
        /// For α, β > 0, skewness = (2(β − α) sqrt(α + β + 1)) / ((α + β + 2) sqrt(αβ)).
        ///
        /// - Note: Returns nil if not finite for the chosen numeric type.
        public var skewness: T? { dyn.skewness }

        /// The kurtosis of the distribution (Pearson’s definition).
        ///
        /// - Note: Returns nil if not finite for the chosen numeric type.
        public var kurtosis: T? { dyn.kurtosis }

        /// The excess kurtosis of the distribution (`kurtosis − 3`).
        ///
        /// - Note: Returns nil if not finite for the chosen numeric type.
        public var kurtosisExcess: T? { dyn.kurtosisExcess }

        // MARK: - Reliability functions

        /// The hazard (failure) rate function h(x) = f(x) / (1 − F(x)).
        ///
        /// - Parameter x: Evaluation point.
        /// - Returns: The hazard at `x`.
        /// - Throws: DistributionError from the underlying implementation if evaluation fails.
        public func hazard(_ x: T) throws -> T { try dyn.hazard(x) }

        /// The cumulative hazard function H(x) = −log(1 − F(x)).
        ///
        /// - Parameter x: Evaluation point.
        /// - Returns: The cumulative hazard at `x`.
        /// - Throws: DistributionError from the underlying implementation if evaluation fails.
        public func chf(_ x: T) throws -> T { try dyn.chf(x) }
        
        // MARK: - Lattice metadata (not applicable)

        /// The lattice step size, if the distribution is lattice-supported.
        ///
        /// Always nil for this continuous distribution.
        public var latticeStep: T? { nil }

        /// The lattice origin, if the distribution is lattice-supported.
        ///
        /// Always nil for this continuous distribution.
        public var latticeOrigin: T? { nil }

        // MARK: - Entropy

        /// The differential entropy of the beta distribution (in nats).
        ///
        /// Formula:
        /// - H = ln B(α, β)
        ///       − (α − 1) ψ(α)
        ///       − (β − 1) ψ(β)
        ///       + (α + β − 2) ψ(α + β)
        ///
        /// where B is the beta function and ψ is the digamma function.
        ///
        /// - Returns: The entropy as `T`, or `nil` if not finite for the chosen numeric type
        ///   or if a failure occurs.
        ///
        /// Implementation details:
        /// - Uses SpecialFunctions.logGamma and SpecialFunctions.digamma to compute log B(α, β)
        ///   and digamma terms.
        public var entropy: T? {
            guard alpha > 0 else { return nil }
            guard beta > 0 else { return nil }
            do {
                let logB1 = try SpecialFunctions.logGamma(self.alpha)
                let logB2 = try SpecialFunctions.logGamma(self.beta)
                let logAb2 = try SpecialFunctions.logGamma(self.alpha + self.beta)
                let logB = logB1 + logB2 - logAb2
                let dA = try SpecialFunctions.digamma(self.alpha)
                let dB = try SpecialFunctions.digamma(self.beta)
                let dAB = try SpecialFunctions.digamma(self.alpha + self.beta)
                let term1 = -(self.alpha - 1.0) * dA
                let term2 = term1 - (self.beta - 1.0) * dB
                let term3 = term2 + (self.alpha + self.beta - 2.0) * dAB
                let h = logB + term3
                return h
            }
            catch {
                return nil
            }
        }
        
        // MARK: - Parameter solvers (helpers)

        /// Computes α given a desired mean and variance of a Beta(α, β) distribution.
        ///
        /// Relationships:
        /// - mean = α / (α + β)
        /// - variance = αβ / ((α + β)^2 (α + β + 1))
        ///
        /// This helper solves for α numerically via the Boost bridge, given mean and variance.
        ///
        /// - Parameters:
        ///   - mean: Target mean in [0, 1].
        ///   - variance: Target variance in [0, 0.25], subject to feasibility for a beta distribution.
        /// - Returns: The inferred α as T.
        ///
        /// Warning:
        /// - No explicit validation is performed here; invalid combinations may yield nonsensical values.
        ///   Consider validating (mean, variance) against beta constraints before use.
        public static func find_alpha( mean: T, variance: T) -> T {
            if T.self == Double.self {
                return T(bs_beta_find_alpha_d(Double(mean), Double(variance)))
            }
            if T.self == Float.self {
                return T(bs_beta_find_alpha_f(Float(mean), Float(variance)))
            }
            #if arch(x86_64) || arch(i386)
            return T(bs_beta_find_alpha_l(Float80(mean), Float80(variance)))
            #else
            return T(bs_beta_find_alpha_d(Double(mean), Double(variance)))
            #endif
        }
        
        /// Computes β given a desired mean and variance of a Beta(α, β) distribution.
        ///
        /// See `find_alpha(mean:variance:)` for the relationships among mean, variance, and parameters.
        ///
        /// - Parameters:
        ///   - mean: Target mean in [0, 1].
        ///   - variance: Target variance in [0, 0.25], subject to feasibility.
        /// - Returns: The inferred β as T.
        ///
        /// Warning:
        /// - No explicit validation is performed here; invalid combinations may yield nonsensical values.
        public static func find_beta( mean: T, variance: T) -> T {
            if T.self == Double.self {
                return T(bs_beta_find_beta_d(Double(mean), Double(variance)))
            }
            if T.self == Float.self {
                return T(bs_beta_find_beta_f(Float(mean), Float(variance)))
            }
            #if arch(x86_64) || arch(i386)
            return T(bs_beta_find_beta_l(Float80(mean), Float80(variance)))
            #else
            return T(bs_beta_find_beta_d(Double(mean), Double(variance)))
            #endif
        }
        
        
        /// Solves for α given a fixed β, a point x ∈ [0, 1], and a target probability P(X ≤ x) = probability.
        ///
        /// This is a convenience wrapper around Boost-based numerical inversion routines.
        ///
        /// - Parameters:
        ///   - beta: The fixed shape parameter β > 0.
        ///   - x: The CDF evaluation point in [0, 1].
        ///   - probability: Target CDF value in [0, 1].
        /// - Returns: The inferred α as T.
        ///
        /// Note:
        /// - Convergence and validity depend on the input combination; no explicit validation is performed here.
        public static func find_alpha_from_beta( beta: T, x: T, probability: T) -> T {
            if T.self == Double.self {
                return T(bs_beta_find_alpha_from_beta_d(Double(beta), Double(x), Double(probability)))
            }
            if T.self == Float.self {
                return T(bs_beta_find_alpha_from_beta_f(Float(beta), Float(x), Float(probability)))
            }
            #if arch(x86_64) || arch(i386)
            return T(bs_beta_find_alpha_from_beta_l(Float80(beta), Float80(x), Float80(probability)))
            #else
            return T(
                bs_beta_find_alpha_from_beta_d(
                    Double(beta),
                    Double(x),
                    Double(probability)
                )
            )
            #endif
        }

        /// Solves for β given a fixed α, a point x ∈ [0, 1], and a target probability P(X ≤ x) = probability.
        ///
        /// This is a convenience wrapper around Boost-based numerical inversion routines.
        ///
        /// - Parameters:
        ///   - beta: Represents α > 0 (naming preserved for API stability; this parameter is α).
        ///   - x: The CDF evaluation point in [0, 1].
        ///   - probability: Target CDF value in [0, 1].
        /// - Returns: The inferred β as T.
        ///
        /// Note:
        /// - Convergence and validity depend on the input combination; no explicit validation is performed here.
        public static func find_beta_from_alpha( beta: T, x: T, probability: T) -> T {
            if T.self == Double.self {
                return T(bs_beta_find_beta_from_alpha_d(Double(beta), Double(x), Double(probability)))
            }
            if T.self == Float.self {
                return T(bs_beta_find_beta_from_alpha_f(Float(beta), Float(x), Float(probability)))
            }
            #if arch(x86_64) || arch(i386)
            return T(bs_beta_find_alpha_from_beta_l(Float80(beta), Float80(x), Float80(probability)))
            #else
            return T(
                bs_beta_find_beta_from_alpha_d(
                    Double(beta),
                    Double(x),
                    Double(probability)
                )
            )
            #endif
        }
    }
}
