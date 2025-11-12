//
//  Created by Volker Thieme
//  Copyright © 2025 Volker Thieme.
//  See the project root LICENSE for licensing information.
//

import SwiftyBoostPrelude

extension Distribution {
    /// Inverse gamma distribution Γ⁻¹(α, β) with shape `α > 0` and scale `β > 0`.
    ///
    /// Definitions
    /// - Support: `x > 0`
    /// - PDF: f(x; α, β) = β^α / Γ(α) · x^(−α−1) · exp(−β / x)
    /// - CDF (lower tail): F(x) = Q(α, β / x)
    /// - SF (upper tail): S(x) = P(α, β / x)
    /// - Quantile: Uses the inverse regularized gamma functions via Boost
    ///
    /// Moments (where defined)
    /// - Mean = β / (α − 1) for α > 1, else undefined
    /// - Variance = β² / [(α − 1)² (α − 2)] for α > 2
    /// - Mode = β / (α + 1)
    /// - Skewness = 4 √(α − 2) / (α − 3) for α > 3
    /// - Excess kurtosis = (30 α − 66) / [(α − 3)(α − 4)] for α > 4
    ///
    /// This wrapper keeps a lightweight handle to the Boost.Math implementation and
    /// delegates all evaluation to ``Distribution/Dynamic`` for consistent behaviour
    /// across floating-point precisions.
    public struct InverseGamma<T: Real & BinaryFloatingPoint & Sendable>: Sendable, DistributionProtocol {
        public typealias RealType = T

        /// Shape parameter α (> 0).
        public let shape: T
        /// Scale parameter β (> 0).
        public let scale: T

        private let dyn: Distribution.Dynamic<T>

        /// Creates an inverse gamma distribution Γ⁻¹(α, β).
        ///
        /// - Parameters:
        ///   - shape: Shape α (> 0).
        ///   - scale: Scale β (> 0). Defaults to 1.
        /// - Throws:
        ///   - ``DistributionError/parameterNotFinite(name:value:)`` if a parameter is non-finite.
        ///   - ``DistributionError/parameterNotPositive(name:value:)`` if a parameter is not strictly positive.
        public init(shape: T, scale: T = 1) throws {
            guard shape.isFinite else {
                throw DistributionError<T>.parameterNotFinite(name: "shape", value: shape)
            }
            guard scale.isFinite else {
                throw DistributionError<T>.parameterNotFinite(name: "scale", value: scale)
            }
            guard shape > 0 else {
                throw DistributionError<T>.parameterNotPositive(name: "shape", value: shape)
            }
            guard scale > 0 else {
                throw DistributionError<T>.parameterNotPositive(name: "scale", value: scale)
            }
            self.shape = shape
            self.scale = scale
            self.dyn = try Distribution.Dynamic<T>(
                distributionName: "inverse_gamma",
                parameters: [
                    "shape": shape,
                    "scale": scale,
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

        // MARK: Quantiles

        public func quantile(_ p: T) throws -> T { try dyn.quantile(p) }
        public func quantileComplement(_ q: T) throws -> T { try dyn.quantileComplement(q) }

        // MARK: Moments & summary statistics

        public var mean: T? { dyn.mean }
        public var variance: T? { dyn.variance }
        public var mode: T? { dyn.mode }
        public var median: T { dyn.median }
        public var skewness: T? { dyn.skewness }
        public var kurtosis: T? { dyn.kurtosis }
        public var kurtosisExcess: T? { dyn.kurtosisExcess }
        public var entropy: T? { dyn.entropy }

        // MARK: Discrete metadata

        public var latticeStep: T? { nil }
        public var latticeOrigin: T? { nil }

        // MARK: KL divergence

        public func klDivergence<D>(
            relativeTo other: D,
            options: Distribution.KLDivergenceOptions<T>
        ) throws -> T? where D: DistributionProtocol, D.RealType == T {
            try DistributionKLDivergenceHelper.evaluate(lhs: self, rhs: other, options: options)
        }
public var isDiscrete: Bool { dyn.isDiscrete }
    }
}
