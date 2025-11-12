//
//  Created by Volker Thieme
//  Copyright © 2025 Volker Thieme.
//  See the project root LICENSE for licensing information.
//

import SwiftyBoostPrelude

extension Distribution {
    /// Inverse normal (inverse Gaussian / Wald) distribution with positive mean `μ`
    /// and shape parameter `λ`.
    ///
    /// Definitions
    /// - Support: `x > 0`
    /// - PDF: f(x; μ, λ) = √(λ / (2πx³)) · exp(−λ(x − μ)² / (2 μ² x))
    /// - Mean = μ
    /// - Variance = μ³ / λ
    /// - Mode = μ (√(1 + 9 μ² / (4 λ²)) − 3 μ / (2 λ))
    /// - Leptokurtic with skewness 3 √(μ / λ) and excess kurtosis 15 μ / λ
    ///
    /// The implementation delegates to the Boost.Math inverse Gaussian distribution
    /// through the dynamic factory, providing unified numerics across Float, Double,
    /// and Float80 (where available).
    public struct InverseNormal<T: Real & BinaryFloatingPoint & Sendable>: Sendable, DistributionProtocol {
        public typealias RealType = T

        /// Distribution mean μ (> 0).
        public let mu: T
        /// Shape parameter λ (> 0).
        public let lambda: T

        private let dyn: Distribution.Dynamic<T>

        /// Creates an inverse normal (inverse Gaussian) distribution.
        ///
        /// - Parameters:
        ///   - mean: Mean μ (> 0). Defaults to 1.
        ///   - shape: Shape parameter λ (> 0). Defaults to 1.
        /// - Throws:
        ///   - ``DistributionError/parameterNotFinite(name:value:)`` if a parameter is non-finite.
        ///   - ``DistributionError/parameterNotPositive(name:value:)`` if a parameter is non-positive.
        public init(mean: T = 1, shape: T = 1) throws {
            guard mean.isFinite else {
                throw DistributionError<T>.parameterNotFinite(name: "mean", value: mean)
            }
            guard shape.isFinite else {
                throw DistributionError<T>.parameterNotFinite(name: "shape", value: shape)
            }
            guard mean > 0 else {
                throw DistributionError<T>.parameterNotPositive(name: "mean", value: mean)
            }
            guard shape > 0 else {
                throw DistributionError<T>.parameterNotPositive(name: "shape", value: shape)
            }
            self.mu = mean
            self.lambda = shape
            self.dyn = try Distribution.Dynamic<T>(
                distributionName: "inverse_gaussian",
                parameters: [
                    "mean": mean,
                    "lambda": shape,
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

    /// Synonym for ``Distribution/InverseNormal`` emphasising the inverse Gaussian terminology.
    public typealias InverseGaussian<T: Real & BinaryFloatingPoint & Sendable> = InverseNormal<T>
    public typealias WaldDistribtion<T: Real & BinaryFloatingPoint & Sendable> = InverseNormal<T>
}
