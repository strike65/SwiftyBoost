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
    /// Fisher–Snedecor F distribution.
    ///
    /// This distribution arises as the ratio of two scaled chi-square random variables
    /// and is commonly used in analysis of variance (ANOVA) and hypothesis testing.
    ///
    /// - Parameters:
    ///   - T: The real scalar type used for computation. Supports `Double`, `Float`,
    ///        and on some architectures `Float80`. The implementation dispatches to
    ///        the most appropriate C backend for the chosen precision.
    ///
    /// The support of the F distribution is `[0, +∞)`. Summary statistics are
    /// only defined for certain ranges of the degrees of freedom; when undefined
    /// or numerically non-finite in the underlying backend, optional properties
    /// such as `mean`, `variance`, `skewness`, and `kurtosis` will return `nil`.
    struct FisherF<T: BinaryFloatingPoint & Sendable>: Sendable, DistributionProtocol {
        typealias Real = T

        // Internal RAII box for the C-allocated distribution object.
        private final class Box: @unchecked Sendable {
            let raw: UnsafeRawPointer
            let free: (UnsafeMutableRawPointer) -> Void
            init(raw: UnsafeRawPointer, free: @escaping (UnsafeMutableRawPointer) -> Void) { self.raw = raw; self.free = free }
            deinit { free(UnsafeMutableRawPointer(mutating: raw)) }
        }

        /// Numerator degrees of freedom (df1). Must be positive and finite.
        public let degreesOfFreedom1: T
        /// Denominator degrees of freedom (df2). Must be positive and finite.
        public let degreesOfFreedom2: T
        private let box: Box

        /// Creates a Fisher–Snedecor F distribution with the given degrees of freedom.
        ///
        /// - Parameters:
        ///   - df1: Numerator degrees of freedom. Must be `> 0` and finite.
        ///   - df2: Denominator degrees of freedom. Must be `> 0` and finite.
        ///
        /// - Throws:
        ///   - `DistributionError.parameterNotPositive` if either degree of freedom is not strictly positive.
        ///   - `DistributionError.parameterNotFinite` if either degree of freedom is not finite.
        ///   - `DistributionError.parameterOutOfRange` if the underlying backend rejects the parameters.
        public init(degreesOfFreedom1 df1: T, degreesOfFreedom2 df2: T) throws {
            guard df1 > 0 else { throw DistributionError.parameterNotPositive(name: "degreesOfFreedom1") }
            guard df2 > 0 else { throw DistributionError.parameterNotPositive(name: "degreesOfFreedom2") }
            guard df1.isFinite else { throw DistributionError.parameterNotFinite(name: "degreesOfFreedom1") }
            guard df2.isFinite else { throw DistributionError.parameterNotFinite(name: "degreesOfFreedom2") }
            self.degreesOfFreedom1 = df1
            self.degreesOfFreedom2 = df2
            if T.self == Double.self {
                guard let h = bs_fisher_f_make(Double(df1), Double(df2)) else {
                    throw DistributionError.generalError(msg: "Distribution not initialized")
                }
                self.box = Box(raw: UnsafeRawPointer(h), free: { bs_fisher_f_free($0) })
            } else if T.self == Float.self {
                guard let h = bs_fisher_f_make_f(Float(df1), Float(df2)) else {
                    throw DistributionError.generalError(msg: "Distribution not initialized")
                }
                self.box = Box(raw: UnsafeRawPointer(h), free: { bs_fisher_f_free_f($0) })
            } else {
                #if arch(x86_64) || arch(i386)
                guard let h = bs_fisher_f_make_l(Float80(df1), Float80(df2)) else {
                    throw DistributionError.generalError(msg: "Distribution not initialized")
                }
                self.box = Box(raw: UnsafeRawPointer(h), free: { bs_fisher_f_free_l($0) })
                #else
                guard let h = bs_fisher_f_make(Double(df1), Double(df2)) else {
                    throw DistributionError.generalError(msg: "Distribution not initialized")
                }
                self.box = Box(raw: UnsafeRawPointer(h), free: { bs_fisher_f_free($0) })
                #endif
            }
        }

        /// The lower bound of the support, inclusive.
        ///
        /// For an F distribution, this is typically `0`.
        public var supportLowerBound: T {
            if T.self == Double.self { let r = bs_fisher_f_range(box.raw); return T(r.lower) }
            if T.self == Float.self  { let r = bs_fisher_f_range_f(box.raw); return T(r.lower) }
            #if arch(x86_64) || arch(i386)
            let r = bs_fisher_f_range_l(box.raw); return T(r.lower)
            #else
            let r = bs_fisher_f_range(box.raw);  return T(r.lower)
            #endif
        }

        /// The upper bound of the support.
        ///
        /// For an F distribution, this is typically `+∞`.
        public var supportUpperBound: T {
            if T.self == Double.self { let r = bs_fisher_f_range(box.raw); return T(r.upper) }
            if T.self == Float.self  { let r = bs_fisher_f_range_f(box.raw); return T(r.upper) }
            #if arch(x86_64) || arch(i386)
            let r = bs_fisher_f_range_l(box.raw); return T(r.upper)
            #else
            let r = bs_fisher_f_range(box.raw);  return T(r.upper)
            #endif
        }

        /// The full support range `(lower, upper)`.
        public var range: (lower: T, upper: T) {
            if T.self == Double.self { let r = bs_fisher_f_range(box.raw);   return (T(r.lower), T(r.upper)) }
            if T.self == Float.self  { let r = bs_fisher_f_range_f(box.raw); return (T(r.lower), T(r.upper)) }
            #if arch(x86_64) || arch(i386)
            let r = bs_fisher_f_range_l(box.raw); return (T(r.lower), T(r.upper))
            #else
            let r = bs_fisher_f_range(box.raw);   return (T(r.lower), T(r.upper))
            #endif
        }

        /// Probability density function (PDF) evaluated at `x`.
        ///
        /// - Parameter x: Point at which to evaluate the density. Values outside the support yield 0.
        /// - Returns: The density value at `x`.
        /// - Throws: Propagates any backend errors encountered during evaluation.
        public func pdf(_ x: T) throws -> T {
            if T.self == Double.self { return T(bs_fisher_f_pdf(box.raw, Double(x))) }
            if T.self == Float.self { return T(bs_fisher_f_pdf_f(box.raw, Float(x))) }
            #if arch(x86_64) || arch(i386)
            return T(bs_fisher_f_pdf_l(box.raw, Float80(x)))
            #else
            return T(bs_fisher_f_pdf(box.raw, Double(x)))
            #endif
        }

        /// Natural logarithm of the PDF at `x`.
        ///
        /// This is computed as `log(pdf(x))`.
        /// - Parameter x: Point at which to evaluate the log-density.
        /// - Returns: The log-density at `x`.
        /// - Throws: See `pdf(_:)`.
        public func logPdf(_ x: T) throws -> T {
            let p = try pdf(x)
            return T(log(Double(p)))
        }

        /// Cumulative distribution function (CDF) evaluated at `x`.
        ///
        /// - Parameter x: Upper integration bound.
        /// - Returns: `P(X ≤ x)`.
        /// - Throws: Propagates any backend errors encountered during evaluation.
        public func cdf(_ x: T) throws -> T {
            if T.self == Double.self { return T(bs_fisher_f_cdf(box.raw, Double(x))) }
            if T.self == Float.self { return T(bs_fisher_f_cdf_f(box.raw, Float(x))) }
            #if arch(x86_64) || arch(i386)
            return T(bs_fisher_f_cdf_l(box.raw, Float80(x)))
            #else
            return T(bs_fisher_f_cdf(box.raw, Double(x)))
            #endif
        }

        /// Survival function (SF), i.e. complementary CDF, evaluated at `x`.
        ///
        /// - Parameter x: Threshold.
        /// - Returns: `P(X > x) = 1 - CDF(x)`.
        /// - Throws: Propagates any backend errors encountered during evaluation.
        public func sf(_ x: T) throws -> T {
            if T.self == Double.self { return T(bs_fisher_f_ccdf(box.raw, Double(x))) }
            if T.self == Float.self { return T(bs_fisher_f_ccdf_f(box.raw, Float(x))) }
            #if arch(x86_64) || arch(i386)
            return T(bs_fisher_f_ccdf_l(box.raw, Float80(x)))
            #else
            return T(bs_fisher_f_ccdf(box.raw, Double(x)))
            #endif
        }

        /// Quantile function (inverse CDF).
        ///
        /// - Parameter p: Probability in `[0, 1]`.
        /// - Returns: Smallest `x` such that `CDF(x) ≥ p`.
        /// - Throws: If `p` is outside `[0, 1]` or the backend fails to converge.
        public func quantile(_ p: T) throws -> T {
            if T.self == Double.self { return T(bs_fisher_f_quantile(box.raw, Double(p))) }
            if T.self == Float.self { return T(bs_fisher_f_quantile_f(box.raw, Float(p))) }
            #if arch(x86_64) || arch(i386)
            return T(bs_fisher_f_quantile_l(box.raw, Float80(p)))
            #else
            return T(bs_fisher_f_quantile(box.raw, Double(p)))
            #endif
        }

        /// Complementary quantile function (inverse survival function).
        ///
        /// - Parameter q: Tail probability in `[0, 1]`.
        /// - Returns: Smallest `x` such that `SF(x) ≤ q`, i.e. `P(X > x) ≤ q`.
        /// - Throws: If `q` is outside `[0, 1]` or the backend fails to converge.
        public func quantileComplement(_ q: T) throws -> T {
            if T.self == Double.self { return T(bs_fisher_f_quantile_complement(box.raw, Double(q))) }
            if T.self == Float.self { return T(bs_fisher_f_quantile_complement_f(box.raw, Float(q))) }
            #if arch(x86_64) || arch(i386)
            return T(bs_fisher_f_quantile_complement_l(box.raw, Float80(q)))
            #else
            return T(bs_fisher_f_quantile_complement(box.raw, Double(q)))
            #endif
        }

        /// Mean of the distribution, when defined.
        ///
        /// - Returns: The mean if finite; otherwise `nil` (e.g. undefined when `df2 ≤ 2`).
        public var mean: T? {
            let m: T = {
                if T.self == Double.self { return T(bs_fisher_f_mean(box.raw)) }
                if T.self == Float.self  { return T(bs_fisher_f_mean_f(box.raw)) }
                #if arch(x86_64) || arch(i386)
                return T(bs_fisher_f_mean_l(box.raw))
                #else
                return T(bs_fisher_f_mean(box.raw))
                #endif
            }()
            return m.isFinite ? m : nil
        }

        /// Variance of the distribution, when defined.
        ///
        /// - Returns: The variance if finite; otherwise `nil` (e.g. undefined when `df2 ≤ 4`).
        public var variance: T? {
            let v: T = {
                if T.self == Double.self { return T(bs_fisher_f_variance(box.raw)) }
                if T.self == Float.self  { return T(bs_fisher_f_variance_f(box.raw)) }
                #if arch(x86_64) || arch(i386)
                return T(bs_fisher_f_variance_l(box.raw))
                #else
                return T(bs_fisher_f_variance(box.raw))
                #endif
            }()
            return v.isFinite ? v : nil
        }

        /// Mode of the distribution, when defined.
        ///
        /// - Note: For the F distribution, the mode exists for `df1 > 2`.
        public var mode: T? {
            if T.self == Double.self { return T(bs_fisher_f_mode(box.raw)) }
            if T.self == Float.self  { return T(bs_fisher_f_mode_f(box.raw)) }
            #if arch(x86_64) || arch(i386)
            return T(bs_fisher_f_mode_l(box.raw))
            #else
            return T(bs_fisher_f_mode(box.raw))
            #endif
        }

        /// Median of the distribution.
        ///
        /// - Note: There is no simple closed form; this value is computed numerically by the backend.
        public var median: T {
            if T.self == Double.self { return T(bs_fisher_f_median(box.raw)) }
            if T.self == Float.self  { return T(bs_fisher_f_median_f(box.raw)) }
            #if arch(x86_64) || arch(i386)
            return T(bs_fisher_f_median_l(box.raw))
            #else
            return T(bs_fisher_f_median(box.raw))
            #endif
        }

        /// Skewness of the distribution, when defined.
        ///
        /// - Returns: The skewness if finite; otherwise `nil` (e.g. undefined when `df2 ≤ 6`).
        public var skewness: T? {
            let s: T = {
                if T.self == Double.self { return T(bs_fisher_f_skewness(box.raw)) }
                if T.self == Float.self  { return T(bs_fisher_f_skewness_f(box.raw)) }
                #if arch(x86_64) || arch(i386)
                return T(bs_fisher_f_skewness_l(box.raw))
                #else
                return T(bs_fisher_f_skewness(box.raw))
                #endif
            }()
            return s.isFinite ? s : nil
        }

        /// Kurtosis of the distribution (including the +3 term), when defined.
        ///
        /// - Returns: The kurtosis if finite; otherwise `nil` (e.g. undefined when `df2 ≤ 8`).
        public var kurtosis: T? {
            let k: T = {
                if T.self == Double.self { return T(bs_fisher_f_kurtosis(box.raw)) }
                if T.self == Float.self  { return T(bs_fisher_f_kurtosis_f(box.raw)) }
                #if arch(x86_64) || arch(i386)
                return T(bs_fisher_f_kurtosis_l(box.raw))
                #else
                return T(bs_fisher_f_kurtosis(box.raw))
                #endif
            }()
            return k.isFinite ? k : nil
        }

        /// Excess kurtosis of the distribution (kurtosis − 3), when defined.
        ///
        /// - Returns: The excess kurtosis if finite; otherwise `nil`.
        public var kurtosisExcess: T? {
            let e: T = {
                if T.self == Double.self { return T(bs_fisher_f_kurtosis_excess(box.raw)) }
                if T.self == Float.self  { return T(bs_fisher_f_kurtosis_excess_f(box.raw)) }
                #if arch(x86_64) || arch(i386)
                return T(bs_fisher_f_kurtosis_excess_l(box.raw))
                #else
                return T(bs_fisher_f_kurtosis_excess(box.raw))
                #endif
            }()
            return e.isFinite ? e : nil
        }

        // MARK: Hazards

        /// Hazard (failure) rate at `x`.
        ///
        /// Defined as `h(x) = pdf(x) / SF(x)`.
        /// - Parameter x: Point at which to evaluate the hazard.
        /// - Returns: The hazard rate at `x`.
        /// - Throws: Propagates any backend errors encountered during evaluation.
        public func hazard(_ x: T) throws -> T {
            if T.self == Double.self { return T(bs_fisher_f_hazard(box.raw, Double(x))) }
            if T.self == Float.self  { return T(bs_fisher_f_hazard_f(box.raw, Float(x))) }
            #if arch(x86_64) || arch(i386)
            return T(bs_fisher_f_hazard_l(box.raw, Float80(x)))
            #else
            return T(bs_fisher_f_hazard(box.raw, Double(x)))
            #endif
        }

        /// Cumulative hazard function (CHF) at `x`.
        ///
        /// Defined as `H(x) = -log(SF(x))`.
        /// - Parameter x: Point at which to evaluate the cumulative hazard.
        /// - Returns: The cumulative hazard at `x`.
        /// - Throws: Propagates any backend errors encountered during evaluation.
        public func chf(_ x: T) throws -> T {
            if T.self == Double.self { return T(bs_fisher_f_chf(box.raw, Double(x))) }
            if T.self == Float.self  { return T(bs_fisher_f_chf_f(box.raw, Float(x))) }
            #if arch(x86_64) || arch(i386)
            return T(bs_fisher_f_chf_l(box.raw, Float80(x)))
            #else
            return T(bs_fisher_f_chf(box.raw, Double(x)))
            #endif
        }
        
        // Lattice/discrete-only properties (continuous ⇒ nil)

        /// For continuous distributions, this is `nil`.
        public var latticeStep: T? { nil }

        /// For continuous distributions, this is `nil`.
        public var latticeOrigin: T? { nil }

        /// Differential entropy of the distribution, when numerically available.
        ///
        /// The entropy is computed via the Beta-parameterization using
        /// `a = df1 / 2`, `b = df2 / 2`, and a log-scale correction by `log(df2/df1)`.
        /// Returns `nil` if parameters are invalid or the special functions fail.
        public var entropy: T? {
            let a: T = self.degreesOfFreedom1 / T(2)
            let b: T = self.degreesOfFreedom2 / T(2)
            let s: T = self.degreesOfFreedom2 / self.degreesOfFreedom1
            guard a > 0, b > 0 else { return nil }
            do {
                let lgA: T = try SpecialFunctions.logGamma(a)
                let lgB: T = try SpecialFunctions.logGamma(b)
                let lgAB: T = try SpecialFunctions.logGamma(a + b)
                let psiA: T = try SpecialFunctions.digamma(a)
                let psiB: T = try SpecialFunctions.digamma(b)
                let psiAB: T = try SpecialFunctions.digamma(a + b)
                let lnB: T =  lgA + lgB - lgAB
                let hBP: T = lnB - (a - 1) * psiA - (b + 1) * psiB + (a + b) * psiAB
                #if arch(x86_64) || arch(i386)
                return hBP + T(log(Float80(s)))
                #else
                return hBP + T(log(Double(s)))
                #endif
            } catch {
                return nil
            }
        }

    }
}
