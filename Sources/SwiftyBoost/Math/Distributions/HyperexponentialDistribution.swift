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
    /// Hyperexponential probability distribution (mixture of exponential phases).
    ///
    /// A hyperexponential distribution models the waiting time for the completion of one
    /// of several exponential phases chosen at random. It is parameterised by:
    ///
    /// - `rates`: The positive rate (λ) for each exponential phase.
    /// - `probabilities`: Optional mixing probabilities for the phases (defaults to a uniform
    ///   mixture when omitted). Probabilities need not be normalised; they are rescaled to sum
    ///   to 1 if supplied.
    ///
    /// Implementation notes:
    /// - This wrapper is backed by the runtime ``Distribution/Dynamic`` factory,
    ///   which bridges into Boost.Math’s `hyperexponential_distribution`.
    /// - Parameters are validated on construction (rates must be strictly positive and finite;
    ///   probabilities must be non-negative with at least one positive entry).
    /// - Hazard-related quantities are computed via the bridged backend.
    public struct Hyperexponential<T: Real & BinaryFloatingPoint & Sendable>: Sendable, DistributionProtocol {
        /// The underlying floating-point type.
        public typealias RealType = T

        /// Mixing probabilities for each exponential phase (normalised to sum to 1).
        public let probabilities: [T]

        /// Rate parameters (λᵢ) for each exponential phase.
        public let rates: [T]

        /// Number of exponential phases in the mixture.
        public var phaseCount: Int { rates.count }

        /// Backing dynamic distribution.
        private let dyn: Distribution.Dynamic<T>

        /// Creates a hyperexponential distribution with optional custom mixing probabilities.
        ///
        /// - Parameters:
        ///   - probabilities: Optional array of non-negative weights for each phase.
        ///     If `nil`, phases are mixed uniformly. If provided, the array must contain
        ///     the same number of elements as `rates`, contain at least one positive entry,
        ///     and will be normalised internally to sum to 1.
        ///   - rates: Array of strictly positive rate parameters λᵢ (one per phase).
        ///
        /// - Throws:
        ///   - ``DistributionError/invalidCombination(message:value:)`` if the phase count is zero.
        ///   - ``DistributionError/parameterNotPositive(name:value:)`` if any rate is non-positive.
        ///   - ``DistributionError/parameterNotFinite(name:value:)`` if any rate or probability
        ///     is `NaN`/infinite.
        ///   - ``DistributionError/invalidCombination(message:value:)`` if provided probabilities
        ///     do not match the number of rates or the total weight is zero.
        public init(probabilities: [T]? = nil, rates: [T]) throws {
            guard !rates.isEmpty else {
                throw DistributionError<T>.invalidCombination(
                    message: "Hyperexponential distribution requires at least one rate",
                    value: nil
                )
            }

            // Validate rates
            for (index, rate) in rates.enumerated() {
                guard rate.isFinite else {
                    throw DistributionError<T>.parameterNotFinite(name: "rate[\(index)]", value: rate)
                }
                guard rate > 0 else {
                    throw DistributionError<T>.parameterNotPositive(name: "rate[\(index)]", value: rate)
                }
            }

            self.rates = rates

            // Prepare probabilities
            let normalisedProbabilities: [T]
            if let probs = probabilities {
                guard probs.count == rates.count else {
                    throw DistributionError<T>.invalidCombination(
                        message: "Probabilities must match the number of rates",
                        value: nil
                    )
                }

                var totalWeight: T = 0
                normalisedProbabilities = try probs.enumerated().map { (index, weight) in
                    guard weight.isFinite else {
                        throw DistributionError<T>.parameterNotFinite(name: "probabilities[\(index)]", value: weight)
                    }
                    guard weight >= 0 else {
                        throw DistributionError<T>.parameterNotPositive(name: "probabilities[\(index)]", value: weight)
                    }
                    totalWeight += weight
                    return weight
                }

                guard totalWeight > 0 else {
                    throw DistributionError<T>.invalidCombination(
                        message: "Probabilities must contain at least one positive entry",
                        value: nil
                    )
                }

                self.probabilities = normalisedProbabilities.map { $0 / totalWeight }
            } else {
                let weight = T(1) / T(rates.count)
                self.probabilities = Array(repeating: weight, count: rates.count)
            }

            // Bridge parameters for the dynamic backend.
            var params: [String: T] = [:]
            params.reserveCapacity(rates.count * 2)
            for (idx, rate) in rates.enumerated() {
                params["rate\(idx)"] = rate
            }
            for (idx, probability) in self.probabilities.enumerated() {
                params["prob\(idx)"] = probability
            }

            self.dyn = try Distribution.Dynamic<T>(
                distributionName: "hyperexponential",
                parameters: params
            )
        }

        // MARK: DistributionProtocol — Support

        /// Lower bound of the support (always `0` for mixtures of exponentials).
        public var supportLowerBound: T { dyn.supportLowerBound }

        /// Upper bound of the support (typically `+∞` for mixtures of exponentials).
        public var supportUpperBound: T { dyn.supportUpperBound }

        /// Convenience tuple containing both support bounds.
        public var range: (lower: T, upper: T) { dyn.range }

        // MARK: PDF/CDF/SF/Quantile

        /// Probability density function evaluated at `x`.
        ///
        /// - Parameter x: Evaluation point in `[0, ∞)`.
        /// - Returns: `Σᵢ pᵢ λᵢ e^(−λᵢ x)`.
        public func pdf(_ x: T) throws -> T { try dyn.pdf(x) }

        /// Natural logarithm of the PDF evaluated at `x`.
        ///
        /// - Parameter x: Evaluation point.
        /// - Returns: `log(pdf(x))`.
        public func logPdf(_ x: T) throws -> T { try dyn.logPdf(x) }

        /// Cumulative distribution function `F(x)`.
        ///
        /// - Parameter x: Evaluation point.
        /// - Returns: Probability `P(X ≤ x)`.
        public func cdf(_ x: T) throws -> T { try dyn.cdf(x) }

        /// Survival function `S(x) = P(X > x)`.
        ///
        /// - Parameter x: Evaluation point.
        public func sf(_ x: T) throws -> T { try dyn.sf(x) }

        /// Lower-tail quantile (inverse CDF).
        ///
        /// - Parameter p: Probability in `[0, 1]`.
        public func quantile(_ p: T) throws -> T { try dyn.quantile(p) }

        /// Upper-tail quantile (inverse survival function).
        ///
        /// - Parameter q: Tail probability in `[0, 1]`.
        public func quantileComplement(_ q: T) throws -> T {
            try dyn.quantileComplement(q)
        }

        // MARK: Moments

        public var mean: T? { dyn.mean }

        public var variance: T? { dyn.variance }

        public var skewness: T? { dyn.skewness }

        public var kurtosis: T? { dyn.kurtosis }

        public var kurtosisExcess: T? { dyn.kurtosisExcess }

        public var mode: T? { dyn.mode }

        public var median: T { dyn.median }

        // MARK: Hazards

        /// Instantaneous hazard function `h(x) = f(x) / S(x)`.
        ///
        /// - Parameter x: Evaluation point.
        public func hazard(_ x: T) throws -> T { try dyn.hazard(x) }

        /// Cumulative hazard function `H(x) = −log S(x)`.
        ///
        /// - Parameter x: Evaluation point.
        public func chf(_ x: T) throws -> T { try dyn.chf(x) }

        // MARK: Discrete-only metadata (continuous ⇒ nil)

        /// Lattice spacing for discrete distributions.
        ///
        /// Hyperexponential mixtures are continuous, so this is always `nil`.
        public var latticeStep: T? { nil }

        /// Lattice origin for discrete distributions.
        ///
        /// Hyperexponential mixtures are continuous, so this is always `nil`.
        public var latticeOrigin: T? { nil }

        // MARK: Entropy and divergence

        /// Differential entropy `h[X]` (nats), when available from the backend.
        ///
        /// - Returns: The entropy supplied by Boost.Math, or `nil` when undefined.
        public var entropy: T? { dyn.entropy }

        /// Indicates whether this distribution is discrete (`true`) or continuous (`false`).
        ///
        /// Hyperexponential mixtures are continuous on `[0, ∞)`, so this always returns `false`.
        public var isDiscrete: Bool { dyn.isDiscrete }

        /// Computes the Kullback–Leibler divergence `D_KL(self || other)` when defined.
        ///
        /// - Parameters:
        ///   - other: The reference hyperexponential distribution *Q*.
        ///   - options: Quadrature configuration; defaults to ``Distribution/KLDivergenceOptions/automatic()``.
        /// - Returns: The divergence in nats, or `nil` if it cannot be evaluated.
        /// - Throws: Rethrows any backend or quadrature errors.
        public func klDivergence<D>(
            relativeTo other: D,
            options: Distribution.KLDivergenceOptions<T>
        ) throws -> T? where D: DistributionProtocol, D.RealType == T {
            try DistributionKLDivergenceHelper.evaluate(lhs: self, rhs: other, options: options)
        }
}
}
