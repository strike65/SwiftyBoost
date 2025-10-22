import SwiftyBoostPrelude

extension Distribution {
    /// The Bernoulli distribution modelling a single binary trial with success probability `p`.
    ///
    /// Definition
    /// - Support: discrete values {0, 1}
    /// - Parameter: `p ∈ [0, 1]` is the probability of observing `1` (success).
    /// - PMF: `P(X = 1) = p`, `P(X = 0) = 1 − p`.
    ///
    /// Notes
    /// - Internally delegates to the dynamic distribution factory backed by Boost.Math.
    /// - Conforms to ``DistributionProtocol`` so it shares the same PDF/PMF, CDF, quantile,
    ///   and summary statistic APIs as the continuous distributions.
    /// - `Float`, `Double`, and (on x86_64) `Float80` are supported.
    public struct Bernoulli<T: Real & BinaryFloatingPoint & Sendable>: Sendable, DistributionProtocol {
        /// The floating-point type used by this distribution.
        typealias RealType = T

        /// Probability of success (`P(X = 1)`), stored as provided during initialization.
        public let success_fraction: T
        private let dyn: Distribution.Dynamic<T>

        /// Creates a Bernoulli distribution with success probability `p`.
        ///
        /// - Parameter p: Probability of success in `[0, 1]`.
        /// - Throws: ``DistributionError/invalidCombination(message:)`` when `p` lies outside `[0, 1]`.
        public init(p: T) throws {
            guard p >= 0, p <= 1 else {
                throw DistributionError.invalidCombination(message: "p must be in [0,1]")
            }
            self.success_fraction = p
            self.dyn = try Distribution.Dynamic<T>(
                distributionName: "bernoulli",
                parameters: [
                    "p": p,
                ]
            )
        }

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

        public var mean: T? { dyn.mean }
        public var variance: T? { dyn.variance }
        public var mode: T? { dyn.mode }
        public var median: T { dyn.median }
        public var skewness: T? { dyn.skewness }
        public var kurtosis: T? { dyn.kurtosis }
        public var kurtosisExcess: T? { dyn.kurtosisExcess }

        public func hazard(_ x: T) throws -> T { try dyn.hazard(x) }
        public func chf(_ x: T) throws -> T { try dyn.chf(x) }

        public var latticeStep: T? { nil }
        public var latticeOrigin: T? { nil }

        public var entropy: T? {
            // Entropy is defined only for p in [0, 1]; at the endpoints it is 0.
            guard success_fraction >= 0.0, success_fraction <= 1.0 else { return nil }
            if self.success_fraction.isZero || self.success_fraction == 1 {
                return 0
            }
            // H = −p*log(p) − (1−p)*log(1−p), computed with log1p for stability.
            let p: T = self.success_fraction
            let q: T = 1 - p
            do {
                // log(p) as log1p(p − 1); valid since p ∈ (0, 1) ⇒ p − 1 ∈ (−1, 0)
                let lp: T = try SpecialFunctions.log1p(p - 1)
                // log(1 − p) as log1p(−p); valid since p ∈ (0, 1) ⇒ −p ∈ (−1, 0)
                let lq: T = try SpecialFunctions.log1p(-p)
                return -(p * lp) - q * lq
            } catch {
                return nil
            }
        }
    }
}
