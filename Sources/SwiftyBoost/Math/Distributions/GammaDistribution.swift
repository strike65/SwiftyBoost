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

import CBoostBridge
import Foundation

extension Distribution {
    /// Gamma distribution Γ(k, θ) with shape k > 0 and scale θ > 0.
    ///
    /// Definitions
    /// - PDF: f(x; k, θ) = x^(k−1) e^(−x/θ) / (Γ(k) θ^k), for x ≥ 0
    /// - CDF (lower tail): F(x) = P(k, x/θ)
    /// - SF (upper tail): S(x) = Q(k, x/θ)
    /// - Quantile: F(x) = p ⇒ x = θ · P⁻¹(k, p)
    ///
    /// Moments (where defined)
    /// - Mean = kθ
    /// - Variance = kθ²
    /// - Mode = (k − 1)θ for k ≥ 1; undefined for k < 1
    /// - Skewness = 2 / √k
    /// - Kurtosis (Pearson) = 3 + 6/k
    /// - Excess kurtosis = 6/k
    ///
    /// Notes
    /// - This type constructs the underlying Boost `gamma_distribution` once and
    ///   keeps a small opaque handle; all evaluations reuse this instance.
    /// - Throws on invalid parameters (k ≤ 0 or θ ≤ 0) or invalid arguments to
    ///   evaluators (e.g., x < 0 for PDF/CDF/SF).
    public struct Gamma<T: BinaryFloatingPoint & Sendable>: Sendable, DistributionProtocol {
        typealias Real = T
        // Backing storage: opaque handle + typed free routine.
        private final class Box: @unchecked Sendable {
            let raw: UnsafeRawPointer
            let free: (UnsafeMutableRawPointer) -> Void
            init(raw: UnsafeRawPointer, free: @escaping (UnsafeMutableRawPointer) -> Void) {
                self.raw = raw
                self.free = free
            }
            deinit { free(UnsafeMutableRawPointer(mutating: raw)) }
        }
        
        /// Shape parameter k (> 0).
        public let shape: T
        /// Scale parameter θ (> 0).
        public let scale: T
        private let box: Box
        
        /// Initialize Γ(k, θ).
        ///
        /// - Parameters:
        ///   - shape: k > 0
        ///   - scale: θ > 0 (default 1)
        /// - Throws: `DistributionError.parameterNotPositive` if parameters are invalid.
        public init(shape: T, scale: T = 1) throws {
            guard shape > 0 else { throw DistributionError.parameterNotPositive(name: "shape") }
            guard scale > 0 else { throw DistributionError.parameterNotPositive(name: "scale") }
            self.shape = shape
            self.scale = scale
            
            if T.self == Double.self {
                guard let h = bs_gamma_make_d(Double(shape), Double(scale)) else {
                    throw DistributionError.generalError(msg: "Distribution not initialized")
                }
                self.box = Box(raw: UnsafeRawPointer(h), free: { bs_gamma_free_d($0) })
            } else if T.self == Float.self {
                guard let h = bs_gamma_make_f(Float(shape), Float(scale)) else {
                    throw DistributionError.generalError(msg: "Distribution not initialized")
                }
                self.box = Box(raw: UnsafeRawPointer(h), free: { bs_gamma_free_f($0) })
            } else {
#if arch(x86_64) || arch(i386)
                // Float80 path uses long double backend
                let k = Float80(shape); let th = Float80(scale)
                guard let h = bs_gamma_make_l(k, th) else {
                    throw DistributionError.generalError(msg: "Distribution not initialized")
                }
                self.box = Box(raw: UnsafeRawPointer(h), free: { bs_gamma_free_l($0) })
#else
                // Fallback via Double handle when Float80 is unavailable
                guard let h = bs_gamma_make_d(Double(shape), Double(scale)) else {
                    throw DistributionError.generalError(msg: "Distribution not initialized")
                }
                self.box = Box(raw: UnsafeRawPointer(h), free: { bs_gamma_free_d($0) })
#endif
            }
        }
        
        // MARK: DistributionProtocol — Support
        
        /// The lower bound of the support. For Γ(k, θ), this is 0.
        public var supportLowerBound: T {
            if T.self == Double.self {
                let r = bs_gamma_support_d(box.raw)
                return T(r.lower)
            }
            if T.self == Float.self {
                let r = bs_gamma_support_f(box.raw)
                return T(r.lower)
            }
            #if arch(x86_64) || arch(i386)
            let r = bs_gamma_support_l(box.raw)
            return T(r.lower)
            #else
            let r = bs_gamma_support_d(box.raw)
            return T(r.lower)
            #endif
        }
        
        /// The upper bound of the support. For Γ(k, θ), this is +∞.
        public var supportUpperBound: T {
            if T.self == Double.self {
                let r = bs_gamma_support_d(box.raw)
                return T(r.upper)
            }
            if T.self == Float.self {
                let r = bs_gamma_support_f(box.raw)
                return T(r.upper)
            }
            #if arch(x86_64) || arch(i386)
            let r = bs_gamma_support_l(box.raw)
            return T(r.upper)
            #else
            let r = bs_gamma_support_d(box.raw)
            return T(r.upper)
            #endif
        }

        /// A convenience tuple of the distribution’s overall range.
        ///
        /// For Γ(k, θ) this is typically (0, +∞), possibly adjusted by Boost’s
        /// internal policies for numerical underflow/overflow reporting.
        public var range: (lower: T, upper: T) {
            if T.self == Double.self {
                let r = bs_gamma_range_d(box.raw)
                return (T(r.lower), T(r.upper))
            }
            if T.self == Float.self {
                let r = bs_gamma_range_f(box.raw)
                return (T(r.lower), T(r.upper))
            }
            #if arch(x86_64) || arch(i386)
            let r = bs_gamma_range_l(box.raw)
            return (T(r.lower), T(r.upper))
            #else
            let r = bs_gamma_range_d(box.raw)
            return (T(r.lower), T(r.upper))
            #endif
        }

        // MARK: PDF/CDF/SF/Quantile
        
        /// Evaluate the probability density function f(x) at `x`.
        ///
        /// - Parameter x: Must satisfy x ≥ 0 for Γ(k, θ).
        /// - Returns: The density value at `x`.
        /// - Throws: `DistributionError.parameterOutOfRange` if `x < 0`.
        public func pdf(_ x: T) throws -> T {
            guard x >= 0 else { throw DistributionError.parameterOutOfRange(name: "x", min: 0.0, max: Double.infinity) }
            if T.self == Double.self { return T(bs_gamma_pdf_d(box.raw, Double(x))) }
            if T.self == Float.self { return T(bs_gamma_pdf_f(box.raw, Float(x))) }
#if arch(x86_64) || arch(i386)
            return T(bs_gamma_pdf_l(box.raw, Float80(x)))
#else
            return T(bs_gamma_pdf_d(box.raw, Double(x)))
#endif
        }
        
        /// Evaluate the lower-tail CDF F(x) = P(X ≤ x) at `x`.
        ///
        /// - Parameter x: Must satisfy x ≥ 0 for Γ(k, θ).
        /// - Returns: The lower-tail probability in [0, 1].
        /// - Throws: `DistributionError.parameterOutOfRange` if `x < 0`.
        public func cdf(_ x: T) throws -> T {
            guard x >= 0 else { throw DistributionError.parameterOutOfRange(name: "x", min: 0.0, max: Double.infinity) }
            if x == 0 { return 0 }
            if T.self == Double.self { return T(bs_gamma_cdf_d(box.raw, Double(x))) }
            if T.self == Float.self { return T(bs_gamma_cdf_f(box.raw, Float(x))) }
#if arch(x86_64) || arch(i386)
            return T(bs_gamma_cdf_l(box.raw, Float80(x)))
#else
            return T(bs_gamma_cdf_d(box.raw, Double(x)))
#endif
        }
        
        /// Evaluate the survival function S(x) = P(X > x) at `x`.
        ///
        /// - Parameter x: Must satisfy x ≥ 0 for Γ(k, θ).
        /// - Returns: The upper-tail probability in [0, 1].
        /// - Throws: `DistributionError.parameterOutOfRange` if `x < 0`.
        public func sf(_ x: T) throws -> T {
            guard x >= 0 else { throw DistributionError.parameterOutOfRange(name: "x", min: 0.0, max: Double.infinity) }
            if x == 0 { return 1 }
            if T.self == Double.self { return T(bs_gamma_ccdf_d(box.raw, Double(x))) }
            if T.self == Float.self { return T(bs_gamma_ccdf_f(box.raw, Float(x))) }
#if arch(x86_64) || arch(i386)
            return T(bs_gamma_ccdf_l(box.raw, Float80(x)))
#else
            return T(bs_gamma_ccdf_d(box.raw, Double(x)))
#endif
        }
        
        /// Quantile (inverse lower-tail CDF): returns x such that F(x) = p.
        ///
        /// - Parameter p: Probability in [0, 1].
        /// - Returns: The lower-tail quantile x.
        /// - Throws: `DistributionError.parameterOutOfRange` if `p ∉ [0, 1]`.
        public func quantile(_ p: T) throws -> T {
            guard p >= 0 && p <= 1 else { throw DistributionError.parameterOutOfRange(name: "p", min: 0.0, max: 1.0) }
            if p == 0 { return 0 }
            if p == 1 { return T.infinity }
            if T.self == Double.self { return T(bs_gamma_quantile_d(box.raw, Double(p))) }
            if T.self == Float.self { return T(bs_gamma_quantile_f(box.raw, Float(p))) }
#if arch(x86_64) || arch(i386)
            return T(bs_gamma_quantile_l(box.raw, Float80(p)))
#else
            return T(bs_gamma_quantile_d(box.raw, Double(p)))
#endif
        }
        
        /// Upper-tail quantile (inverse SF): returns x such that S(x) = q.
        ///
        /// - Parameter q: Upper-tail probability in [0, 1].
        /// - Returns: The upper-tail quantile x.
        /// - Throws: `DistributionError.parameterOutOfRange` if `q ∉ [0, 1]`.
        public func quantileComplement(_ q: T) throws -> T {
            guard q >= 0 && q <= 1 else { throw DistributionError.parameterOutOfRange(name: "q", min: 0.0, max: 1.0) }
            if q == 1 { return 0 }
            if q == 0 { return T.infinity }
            if T.self == Double.self { return T(bs_gamma_quantile_complement_d(box.raw, Double(q))) }
            if T.self == Float.self { return T(bs_gamma_quantile_complement_f(box.raw, Float(q))) }
#if arch(x86_64) || arch(i386)
            return T(bs_gamma_quantile_complement_l(box.raw, Float80(q)))
#else
            return T(bs_gamma_quantile_complement_d(box.raw, Double(q)))
#endif
        }
        
        // MARK: Moments
        
        /// The mean E[X] = kθ.
        public var mean: T? {
            if T.self == Double.self { return T(bs_gamma_mean_d(box.raw)) }
            if T.self == Float.self { return T(bs_gamma_mean_f(box.raw)) }
#if arch(x86_64) || arch(i386)
            return T(bs_gamma_mean_l(box.raw))
#else
            return T(bs_gamma_mean_d(box.raw))
#endif
        }
        
        /// The variance Var[X] = kθ².
        public var variance: T? {
            if T.self == Double.self { return T(bs_gamma_variance_d(box.raw)) }
            if T.self == Float.self { return T(bs_gamma_variance_f(box.raw)) }
#if arch(x86_64) || arch(i386)
            return T(bs_gamma_variance_l(box.raw))
#else
            return T(bs_gamma_variance_d(box.raw))
#endif
        }
        
        /// The mode, if defined.
        ///
        /// For Γ(k, θ), the mode is (k − 1)θ when k ≥ 1; undefined otherwise.
        /// Returns `nil` when not finite or undefined.
        public var mode: T? {
            let m: T = {
                if T.self == Double.self { return T(bs_gamma_mode_d(box.raw)) }
                if T.self == Float.self { return T(bs_gamma_mode_f(box.raw)) }
#if arch(x86_64) || arch(i386)
                return T(bs_gamma_mode_l(box.raw))
#else
                return T(bs_gamma_mode_d(box.raw))
#endif
            }()
            return m.isFinite ? m : nil
        }
        
        /// The hazard function h(x) = f(x) / S(x) where defined.
        ///
        /// - Parameter x: Must satisfy x ≥ 0 for Γ(k, θ).
        /// - Returns: The instantaneous hazard at `x`.
        /// - Throws: `DistributionError.parameterOutOfRange` if `x < 0`.
        public func hazard(_ x: T) throws -> T {
            guard x >= 0 else { throw DistributionError.parameterOutOfRange(name: "x", min: 0.0, max: Double.infinity) }
            if T.self == Double.self { return T(bs_gamma_hazard_d(box.raw, Double(x))) }
            if T.self == Float.self { return T(bs_gamma_hazard_f(box.raw, Float(x))) }
#if arch(x86_64) || arch(i386)
            return T(bs_gamma_hazard_l(box.raw, Float80(x)))
#else
            return T(bs_gamma_hazard_d(box.raw, Double(x)))
#endif
        }

        /// The cumulative hazard H(x) = −log S(x) where defined.
        ///
        /// - Parameter x: Must satisfy x ≥ 0 for Γ(k, θ).
        /// - Returns: The cumulative hazard at `x`.
        /// - Throws: `DistributionError.parameterOutOfRange` if `x < 0`.
        public func chf(_ x: T) throws -> T {
            guard x >= 0 else { throw DistributionError.parameterOutOfRange(name: "x", min: 0.0, max: Double.infinity) }
            if T.self == Double.self { return T(bs_gamma_chf_d(box.raw, Double(x))) }
            if T.self == Float.self { return T(bs_gamma_chf_f(box.raw, Float(x))) }
#if arch(x86_64) || arch(i386)
            return T(bs_gamma_chf_l(box.raw, Float80(x)))
#else
            return T(bs_gamma_chf_d(box.raw, Double(x)))
#endif
        }

        /// The natural logarithm of the density, log f(x).
        ///
        /// - Parameter x: Must satisfy x ≥ 0 for Γ(k, θ).
        /// - Returns: The log-PDF at `x`.
        /// - Throws: `DistributionError.parameterOutOfRange` if `x < 0`.
        public func logPdf(_ x: T) throws -> T {
            let p = try pdf(x)
            return T(log(Double(p)))
        }

        // MARK: DistributionProtocol — Additional moments
        
        /// The median of the distribution (numerically evaluated by Boost).
        public var median: T {
            if T.self == Double.self { return T(bs_gamma_median_d(box.raw)) }
            if T.self == Float.self { return T(bs_gamma_median_f(box.raw)) }
#if arch(x86_64) || arch(i386)
            return T(bs_gamma_median_l(box.raw))
#else
            return T(bs_gamma_median_d(box.raw))
#endif
        }
        
        /// The skewness γ₁ = E[((X − μ)/σ)^3], where defined.
        ///
        /// For Γ(k, θ), γ₁ = 2 / √k.
        public var skewness: T? {
            if T.self == Double.self { return T(bs_gamma_skewness_d(box.raw)) }
            if T.self == Float.self { return T(bs_gamma_skewness_f(box.raw)) }
#if arch(x86_64) || arch(i386)
            return T(bs_gamma_skewness_l(box.raw))
#else
            return T(bs_gamma_skewness_d(box.raw))
#endif
        }
        
        /// The kurtosis (Pearson’s β₂) = E[((X − μ)/σ)^4], where defined.
        ///
        /// For Γ(k, θ), β₂ = 3 + 6/k. Equals 3 for a normal distribution.
        public var kurtosis: T? {
            if T.self == Double.self { return T(bs_gamma_kurtosis_d(box.raw)) }
            if T.self == Float.self { return T(bs_gamma_kurtosis_f(box.raw)) }
#if arch(x86_64) || arch(i386)
            return T(bs_gamma_kurtosis_l(box.raw))
#else
            return T(bs_gamma_kurtosis_d(box.raw))
#endif
        }
        
        /// The excess kurtosis γ₂ = kurtosis − 3, where defined.
        ///
        /// For Γ(k, θ), γ₂ = 6/k. Equals 0 for a normal distribution.
        public var kurtosisExcess: T? {
            if T.self == Double.self { return T(bs_gamma_kurtosis_excess_d(box.raw)) }
            if T.self == Float.self { return T(bs_gamma_kurtosis_excess_f(box.raw)) }
#if arch(x86_64) || arch(i386)
            return T(bs_gamma_kurtosis_excess_l(box.raw))
#else
            return T(bs_gamma_kurtosis_excess_d(box.raw))
#endif
        }
        
        /// Lattice step for discrete distributions (nil for Γ).
        public var latticeStep: T? { nil }
        /// Lattice origin for discrete distributions (nil for Γ).
        public var latticeOrigin: T? { nil }
        
        /// Entropy H[X] (Shannon), where defined.
        public var entropy: T? {
            if T.self == Double.self { return T(bs_gamma_entropy_d(box.raw)) }
            if T.self == Float.self { return T(bs_gamma_entropy_f(box.raw)) }
#if arch(x86_64) || arch(i386)
            return T(bs_gamma_entropy_l(box.raw))
#else
            return T(bs_gamma_entropy_d(box.raw))
#endif
        }
    }
}
