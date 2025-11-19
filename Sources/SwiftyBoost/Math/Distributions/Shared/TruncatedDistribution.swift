//
//  Created by VT on 19.11.25.
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

/// Errors specific to the truncated distribution wrapper.
public enum TruncationError: Error {
    case valueOutsideTruncationInterval
    case probabilityOutOfBounds
    case invalidNormalizationConstant
}

extension Distribution {
    /// Generic truncated (continuous) distribution wrapping an existing `DistributionProtocol`
    /// and restricting it to an effective interval [supportLowerBound, supportUpperBound].
    ///
    /// The truncation interval is intersected with the base distribution's support.
    /// The base distribution must be continuous (`isDiscrete == false`).
    public struct TruncatedDistribution<Base: DistributionProtocol>: DistributionProtocol {
        public typealias RealType = Base.RealType
        
        // MARK: - Stored properties
        
        /// The underlying base distribution.
        public let base: Base
        
        /// Effective lower bound of the truncated support.
        public let supportLowerBound: RealType
        
        /// Effective upper bound of the truncated support.
        public let supportUpperBound: RealType
        
        /// Precomputed CDF at the lower bound of the truncated interval.
        @usableFromInline
        internal let cdfLower: RealType
        
        /// Precomputed CDF at the upper bound of the truncated interval.
        @usableFromInline
        internal let cdfUpper: RealType
        
        /// Precomputed survival at the lower bound of the truncated interval.
        @usableFromInline
        internal let sfLower: RealType
        
        /// Precomputed survival at the upper bound of the truncated interval.
        @usableFromInline
        internal let sfUpper: RealType
        
        /// Normalization constant Z = F(U) - F(L) = P(L <= X <= U) under the base distribution.
        @usableFromInline
        internal let normalization: RealType
        
        /// log(Z) for numerically stable log-density computations.
        @usableFromInline
        internal let logNormalization: RealType
        
        /// Median of the truncated distribution (precomputed).
        public let median: RealType
        
        // MARK: - Initialization
        
        /// Creates a truncated version of `base` on the interval `[lower, upper]`, intersected
        /// with the base support. At least one of `lower` or `upper` must be non-nil.
        ///
        /// - Parameters:
        ///   - base: The underlying base distribution.
        ///   - lower: Optional lower truncation bound (in original scale).
        ///   - upper: Optional upper truncation bound (in original scale).
        public init(base: Base, lower: RealType? = nil, upper: RealType? = nil) throws {
            precondition(!base.isDiscrete, "TruncatedDistribution currently supports continuous base distributions only.")
            precondition(lower != nil || upper != nil, "At least one truncation bound must be provided.")
            
            let baseLower = base.supportLowerBound
            let baseUpper = base.supportUpperBound
            
            var effectiveLower = baseLower
            var effectiveUpper = baseUpper
            
            if let lowerBound = lower {
                effectiveLower = Swift.max(lowerBound, baseLower)
            }
            if let upperBound = upper {
                effectiveUpper = Swift.min(upperBound, baseUpper)
            }
            
            precondition(effectiveLower < effectiveUpper, "Empty or invalid truncation interval.")
            
            self.base = base
            self.supportLowerBound = effectiveLower
            self.supportUpperBound = effectiveUpper
            
            // Precompute CDF at bounds
            let F_L = try base.cdf(effectiveLower)
            let F_U = try base.cdf(effectiveUpper)
            
            let Z = F_U - F_L
            guard Z > .zero, Z.isFinite else {
                throw TruncationError.invalidNormalizationConstant
            }
            
            self.cdfLower = F_L
            self.cdfUpper = F_U
            self.normalization = Z
            self.logNormalization = RealType.log(Z)
            
            // Precompute survival at bounds (useful for sf / quantileComplement)
            let S_L = try base.sf(effectiveLower)
            let S_U = try base.sf(effectiveUpper)
            self.sfLower = S_L
            self.sfUpper = S_U
            
            // Precompute median as truncated quantile at p = 0.5
            let medianTarget = F_L + RealType(0.5) * Z
            self.median = try base.quantile(medianTarget)
        }
        
        // MARK: - Summary statistics and meta-data
        
        public var range: (lower: RealType, upper: RealType) {
            (lower: supportLowerBound, upper: supportUpperBound)
        }
        
        public var isDiscrete: Bool { false }
        public var latticeStep: RealType? { nil }
        public var latticeOrigin: RealType? { nil }
        
        public var mean: RealType? { nil }
        public var mode: RealType? { nil }
        public var variance: RealType? { nil }
        public var skewness: RealType? { nil }
        public var kurtosis: RealType? { nil }
        public var kurtosisExcess: RealType? { nil }
        public var entropy: RealType? { nil }
        
        // MARK: - Core functions
        
        public func pdf(_ x: RealType) throws -> RealType {
            // Outside truncated support: density is exactly zero.
            if x < supportLowerBound || x > supportUpperBound {
                return .zero
            }
            let baseValue = try base.pdf(x)
            if baseValue <= .zero || !baseValue.isFinite {
                return .zero
            }
            return baseValue / normalization
        }
        
        public func logPdf(_ x: RealType) throws -> RealType {
            // Outside truncated support: log-density is -infinity.
            if x < supportLowerBound || x > supportUpperBound {
                return -.infinity
            }
            let baseLog = try base.logPdf(x)
            return baseLog - logNormalization
        }
        
        public func cdf(_ x: RealType) throws -> RealType {
            // F_tr(x) = 0 for x <= L
            if x <= supportLowerBound {
                return .zero
            }
            // F_tr(x) = 1 for x >= U
            if x >= supportUpperBound {
                return RealType(1)
            }
            // F_tr(x) = (F(x) - F(L)) / Z for L < x < U
            let Fx = try base.cdf(x)
            var value = (Fx - cdfLower) / normalization
            // Guard against small numerical drift
            if value <= .zero { value = .zero }
            if value >= RealType(1) { value = RealType(1) }
            return value
        }
        
        public func sf(_ x: RealType) throws -> RealType {
            // S_tr(x) = 1 for x < L
            if x <= supportLowerBound {
                return RealType(1)
            }
            // S_tr(x) = 0 for x >= U
            if x >= supportUpperBound {
                return .zero
            }
            // S_tr(x) = (S(x) - S(U)) / Z for L <= x < U
            let Sx = try base.sf(x)
            var value = (Sx - sfUpper) / normalization
            if value <= .zero { value = .zero }
            if value >= RealType(1) { value = RealType(1) }
            return value
        }
        
        public func hazard(_ x: RealType) throws -> RealType {
            // Hazard is only defined inside the truncated support.
            if x < supportLowerBound || x > supportUpperBound {
                throw TruncationError.valueOutsideTruncationInterval
            }
            let f = try pdf(x)
            let s = try sf(x)
            if s <= .zero {
                return RealType.infinity
            }
            return f / s
        }
        
        public func chf(_ x: RealType) throws -> RealType {
            // Cumulative hazard is only defined inside the truncated support.
            if x < supportLowerBound || x > supportUpperBound {
                throw TruncationError.valueOutsideTruncationInterval
            }
            let s = try sf(x)
            if s <= .zero {
                return RealType.infinity
            }
            return -RealType.log(s)
        }
        
        // MARK: - Inverses
        
        public func quantile(_ p: RealType) throws -> RealType {
            guard p >= .zero, p <= RealType(1) else {
                throw TruncationError.probabilityOutOfBounds
            }
            if p == .zero {
                return supportLowerBound
            }
            if p == RealType(1) {
                return supportUpperBound
            }
            // F_tr(x) = p  =>  F(x) = F(L) + p * Z
            let target = cdfLower + p * normalization
            return try base.quantile(target)
        }
        
        public func quantileComplement(_ p: RealType) throws -> RealType {
            guard p >= .zero, p <= RealType(1) else {
                throw TruncationError.probabilityOutOfBounds
            }
            if p == .zero {
                return supportUpperBound
            }
            if p == RealType(1) {
                return supportLowerBound
            }
            // S_tr(x) = p  =>  S(x) = S(U) + p * Z
            // Then x solves S(x) = target via the base's upper-tail quantile.
            let target = sfUpper + p * normalization
            return try base.quantileComplement(target)
        }
    }
}
