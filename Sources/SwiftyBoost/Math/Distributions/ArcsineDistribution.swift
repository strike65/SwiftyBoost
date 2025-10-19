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

public extension Distribution {
    /// A continuous arcsine distribution on a finite interval.
    ///
    /// This distribution is supported on the closed interval `[minX, maxX]` and has
    /// the probability density function (PDF)
    ///
    ///     f(x) = 1 / (π * sqrt((x - minX) * (maxX - x)))   for x ∈ (minX, maxX)
    ///
    /// The distribution is U-shaped, with singularities at the endpoints of the support.
    ///
    /// - Note: This type is generic over a floating-point type `T` and supports `Float`, `Double`,
    ///   and, on x86 architectures, `Float80`. Internally it bridges to Boost via `CBoostBridge`.
    struct Arcsine<T: BinaryFloatingPoint & Sendable>: Sendable, DistributionProtocol {
        /// The floating-point type used by this distribution.
        typealias Real = T

        // Internal box to manage the lifetime of the underlying C/Boost object.
        private final class Box: @unchecked Sendable {
            let raw: UnsafeRawPointer
            let free: (UnsafeMutableRawPointer) -> Void
            init(raw: UnsafeRawPointer, free: @escaping (UnsafeMutableRawPointer) -> Void) { self.raw = raw; self.free = free }
            deinit { free(UnsafeMutableRawPointer(mutating: raw)) }
        }

        /// The lower bound of the support (a).
        public let minX: T
        /// The upper bound of the support (b).
        public let maxX: T
        private let box: Box

        /// Creates an arcsine distribution with the given support.
        ///
        /// - Parameters:
        ///   - min_x: The lower bound `a` of the support. Must be finite and strictly less than `max_x`.
        ///   - max_x: The upper bound `b` of the support. Must be finite and strictly greater than `min_x`.
        ///
        /// - Throws: `DistributionError.invalidCombination` if `min_x >= max_x`,
        ///           `DistributionError.parameterNotFinite` if either bound is not finite,
        ///           or `DistributionError.generalError` if the underlying distribution could not be initialized.
        public init(minX min_x: T = 0, maxX max_x: T = 0) throws {
            guard min_x < max_x else { throw DistributionError.invalidCombination(message: "minX must be less than maxX") }
            guard min_x.isFinite else { throw DistributionError.parameterNotFinite(name: "minX") }
            guard max_x.isFinite else { throw DistributionError.parameterNotFinite(name: "maxX") }
            self.minX = min_x
            self.maxX = max_x
            if T.self == Double.self {
                guard let h = bs_arcsine_make(Double(min_x), Double(max_x)) else {
                    throw DistributionError.generalError(msg: "Distribution not initialized")
                }
                self.box = Box(raw: UnsafeRawPointer(h), free: { bs_arcsine_free($0) })
            } else if T.self == Float.self {
                guard let h = bs_arcsine_make_f(Float(min_x), Float(max_x)) else {
                    throw DistributionError.generalError(msg: "Distribution not initialized")
                }
                self.box = Box(raw: UnsafeRawPointer(h), free: { bs_arcsine_free_f($0) })
            } else {
                #if arch(x86_64) || arch(i386)
                guard let h = bs_arcsine_make_l(Float80(df1), Float80(df2)) else {
                    throw DistributionError.generalError(msg: "Distribution not initialized")
                }
                self.box = Box(raw: UnsafeRawPointer(h), free: { bs_arcsine_free_l($0) })
                #else
                guard let h = bs_arcsine_make(Double(min_x), Double(max_x)) else {
                    throw DistributionError.generalError(msg: "Distribution not initialized")
                }
                self.box = Box(raw: UnsafeRawPointer(h), free: { bs_arcsine_free($0) })
                #endif
            }
        }

        /// The lower bound of the support.
        public var supportLowerBound: T {
            if T.self == Double.self { let r = bs_arcsine_range(box.raw); return T(r.lower) }
            if T.self == Float.self  { let r = bs_arcsine_range_f(box.raw); return T(r.lower) }
            #if arch(x86_64) || arch(i386)
            let r = bs_arcsine_range_l(box.raw); return T(r.lower)
            #else
            let r = bs_arcsine_range(box.raw);  return T(r.lower)
            #endif
        }

        /// The upper bound of the support.
        public var supportUpperBound: T {
            if T.self == Double.self { let r = bs_arcsine_range(box.raw); return T(r.upper) }
            if T.self == Float.self  { let r = bs_arcsine_range_f(box.raw); return T(r.upper) }
            #if arch(x86_64) || arch(i386)
            let r = bs_arcsine_range_l(box.raw); return T(r.upper)
            #else
            let r = bs_arcsine_range(box.raw);  return T(r.upper)
            #endif
        }

        /// The support as a (lower, upper) tuple.
        public var range: (lower: T, upper: T) {
            if T.self == Double.self { let r = bs_arcsine_range(box.raw);   return (T(r.lower), T(r.upper)) }
            if T.self == Float.self  { let r = bs_arcsine_range_f(box.raw); return (T(r.lower), T(r.upper)) }
            #if arch(x86_64) || arch(i386)
            let r = bs_arcsine_range_l(box.raw); return (T(r.lower), T(r.upper))
            #else
            let r = bs_arcsine_range(box.raw);   return (T(r.lower), T(r.upper))
            #endif
        }

        /// The probability density function (PDF).
        ///
        /// - Parameter x: Evaluation point.
        /// - Returns: The density at `x`. Values outside the support yield zero.
        public func pdf(_ x: T) throws -> T {
            if T.self == Double.self { return T(bs_arcsine_pdf(box.raw, Double(x))) }
            if T.self == Float.self { return T(bs_arcsine_pdf_f(box.raw, Float(x))) }
            #if arch(x86_64) || arch(i386)
            return T(bs_arcsine_pdf_l(box.raw, Float80(x)))
            #else
            return T(bs_arcsine_pdf(box.raw, Double(x)))
            #endif
        }

        /// The natural logarithm of the PDF.
        ///
        /// - Parameter x: Evaluation point.
        /// - Returns: `log(pdf(x))`. May be `-infinity` at the boundaries.
        public func logPdf(_ x: T) throws -> T {
            let p = try pdf(x)
            return T(log(Double(p)))
        }

        /// The cumulative distribution function (CDF).
        ///
        /// - Parameter x: Evaluation point.
        /// - Returns: `P(X ≤ x)`.
        public func cdf(_ x: T) throws -> T {
            if T.self == Double.self { return T(bs_arcsine_cdf(box.raw, Double(x))) }
            if T.self == Float.self { return T(bs_arcsine_cdf_f(box.raw, Float(x))) }
            #if arch(x86_64) || arch(i386)
            return T(bs_arcsine_cdf_l(box.raw, Float80(x)))
            #else
            return T(bs_arcsine_cdf(box.raw, Double(x)))
            #endif
        }

        /// The survival function (SF), i.e. the complementary CDF.
        ///
        /// - Parameter x: Evaluation point.
        /// - Returns: `P(X > x) = 1 - CDF(x)`.
        public func sf(_ x: T) throws -> T {
            if T.self == Double.self { return T(bs_arcsine_ccdf(box.raw, Double(x))) }
            if T.self == Float.self { return T(bs_arcsine_ccdf_f(box.raw, Float(x))) }
            #if arch(x86_64) || arch(i386)
            return T(bs_arcsine_ccdf_l(box.raw, Float80(x)))
            #else
            return T(bs_arcsine_ccdf(box.raw, Double(x)))
            #endif
        }

        /// The quantile function (inverse CDF).
        ///
        /// - Parameter p: A probability in `[0, 1]`.
        /// - Returns: `x` such that `P(X ≤ x) = p`.
        public func quantile(_ p: T) throws -> T {
            if T.self == Double.self { return T(bs_arcsine_quantile(box.raw, Double(p))) }
            if T.self == Float.self { return T(bs_arcsine_quantile_f(box.raw, Float(p))) }
            #if arch(x86_64) || arch(i386)
            return T(bs_arcsine_quantile_l(box.raw, Float80(p)))
            #else
            return T(bs_arcsine_quantile(box.raw, Double(p)))
            #endif
        }

        /// The complementary quantile function (inverse survival function).
        ///
        /// - Parameter q: A probability in `[0, 1]`.
        /// - Returns: `x` such that `P(X > x) = q`.
        public func quantileComplement(_ q: T) throws -> T {
            if T.self == Double.self { return T(bs_arcsine_quantile_complement(box.raw, Double(q))) }
            if T.self == Float.self { return T(bs_arcsine_quantile_complement_f(box.raw, Float(q))) }
            #if arch(x86_64) || arch(i386)
            return T(bs_arcsine_quantile_complement_l(box.raw, Float80(q)))
            #else
            return T(bs_arcsine_quantile_complement(box.raw, Double(q)))
            #endif
        }

        /// The mean of the distribution.
        ///
        /// - Note: Returns `nil` if not finite for the chosen numeric type.
        public var mean: T? {
            let m: T = {
                if T.self == Double.self { return T(bs_arcsine_mean(box.raw)) }
                if T.self == Float.self  { return T(bs_arcsine_mean_f(box.raw)) }
                #if arch(x86_64) || arch(i386)
                return T(bs_arcsine_mean_l(box.raw))
                #else
                return T(bs_arcsine_mean(box.raw))
                #endif
            }()
            return m.isFinite ? m : nil
        }

        /// The variance of the distribution.
        ///
        /// - Note: Returns `nil` if not finite for the chosen numeric type.
        public var variance: T? {
            let v: T = {
                if T.self == Double.self { return T(bs_arcsine_variance(box.raw)) }
                if T.self == Float.self  { return T(bs_arcsine_variance_f(box.raw)) }
                #if arch(x86_64) || arch(i386)
                return T(bs_arcsine_variance_l(box.raw))
                #else
                return T(bs_arcsine_variance(box.raw))
                #endif
            }()
            return v.isFinite ? v : nil
        }

        /// The mode of the distribution.
        ///
        /// - Note: For the arcsine distribution on a finite interval, the density is unbounded at the endpoints.
        ///   This value reflects the underlying implementation choice in Boost.
        public var mode: T? {
            if T.self == Double.self { return T(bs_arcsine_mode(box.raw)) }
            if T.self == Float.self  { return T(bs_arcsine_mode_f(box.raw)) }
            #if arch(x86_64) || arch(i386)
            return T(bs_arcsine_mode_l(box.raw))
            #else
            return T(bs_arcsine_mode(box.raw))
            #endif
        }

        /// The median of the distribution.
        public var median: T {
            if T.self == Double.self { return T(bs_arcsine_median(box.raw)) }
            if T.self == Float.self  { return T(bs_arcsine_median_f(box.raw)) }
            #if arch(x86_64) || arch(i386)
            return T(bs_arcsine_median_l(box.raw))
            #else
            return T(bs_arcsine_median(box.raw))
            #endif
        }

        /// The skewness of the distribution.
        ///
        /// - Note: Returns `nil` if not finite for the chosen numeric type.
        public var skewness: T? {
            let s: T = {
                if T.self == Double.self { return T(bs_arcsine_skewness(box.raw)) }
                if T.self == Float.self  { return T(bs_arcsine_skewness_f(box.raw)) }
                #if arch(x86_64) || arch(i386)
                return T(bs_arcsine_skewness_l(box.raw))
                #else
                return T(bs_arcsine_skewness(box.raw))
                #endif
            }()
            return s.isFinite ? s : nil
        }

        /// The kurtosis of the distribution.
        ///
        /// - Note: Returns `nil` if not finite for the chosen numeric type.
        public var kurtosis: T? {
            let k: T = {
                if T.self == Double.self { return T(bs_arcsine_kurtosis(box.raw)) }
                if T.self == Float.self  { return T(bs_arcsine_kurtosis_f(box.raw)) }
                #if arch(x86_64) || arch(i386)
                return T(bs_arcsine_kurtosis_l(box.raw))
                #else
                return T(bs_arcsine_kurtosis(box.raw))
                #endif
            }()
            return k.isFinite ? k : nil
        }

        /// The excess kurtosis of the distribution (`kurtosis - 3`).
        ///
        /// - Note: Returns `nil` if not finite for the chosen numeric type.
        public var kurtosisExcess: T? {
            let e: T = {
                if T.self == Double.self { return T(bs_arcsine_kurtosis_excess(box.raw)) }
                if T.self == Float.self  { return T(bs_arcsine_kurtosis_excess_f(box.raw)) }
                #if arch(x86_64) || arch(i386)
                return T(bs_arcsine_kurtosis_excess_l(box.raw))
                #else
                return T(bs_arcsine_kurtosis_excess(box.raw))
                #endif
            }()
            return e.isFinite ? e : nil
        }

        /// The hazard (failure) rate function `h(x) = f(x) / (1 - F(x))`.
        ///
        /// - Parameter x: Evaluation point.
        /// - Returns: The hazard at `x`.
        public func hazard(_ x: T) throws -> T {
            if T.self == Double.self { return T(bs_arcsine_hazard(box.raw, Double(x))) }
            if T.self == Float.self  { return T(bs_arcsine_hazard_f(box.raw, Float(x))) }
            #if arch(x86_64) || arch(i386)
            return T(bs_arcsine_hazard_l(box.raw, Float80(x)))
            #else
            return T(bs_arcsine_hazard(box.raw, Double(x)))
            #endif
        }

        /// The cumulative hazard function `H(x) = -log(1 - F(x))`.
        ///
        /// - Parameter x: Evaluation point.
        /// - Returns: The cumulative hazard at `x`.
        public func chf(_ x: T) throws -> T {
            if T.self == Double.self { return T(bs_arcsine_chf(box.raw, Double(x))) }
            if T.self == Float.self  { return T(bs_arcsine_chf_f(box.raw, Float(x))) }
            #if arch(x86_64) || arch(i386)
            return T(bs_arcsine_chf_l(box.raw, Float80(x)))
            #else
            return T(bs_arcsine_chf(box.raw, Double(x)))
            #endif
        }
        
        /// The lattice step size, if the distribution is lattice-supported.
        ///
        /// Always `nil` for this continuous distribution.
        public var latticeStep: T? { nil }

        /// The lattice origin, if the distribution is lattice-supported.
        ///
        /// Always `nil` for this continuous distribution.
        public var latticeOrigin: T? { nil }

        /// The differential entropy of the arcsine distribution.
        ///
        /// For support width `w = maxX - minX`, the entropy is:
        ///
        ///     H = ln(π / 4) + ln(w) = ln(π w / 4)
        ///
        /// - Returns: The entropy as `T`, or `nil` if not finite for the chosen numeric type.
        public var entropy: T? {
            let width = self.maxX - self.minX
            #if arch(x86_64) || arch(i386)
            let pi: Float80 = Constants.pi
            let qpi: Float80 = Constants.quarterPi
            let res: T = T(qpi) + T(log(Float80(width)))
            return res
            #else
            let pi: Double = Constants.pi
            let qpi: Double = Constants.quarterPi
            let res: T = T(log(qpi)) + T(log(Double(width)))
            return res
            #endif
        }

    }
}
