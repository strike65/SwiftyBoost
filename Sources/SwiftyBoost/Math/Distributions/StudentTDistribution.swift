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
    /// Student's t distribution with ν degrees of freedom (ν > 0).
    ///
    /// Definitions
    /// - PDF: f(x; ν) = Γ((ν+1)/2) / (√(νπ) Γ(ν/2)) · (1 + x²/ν)^{−(ν+1)/2}
    /// - CDF: integral of PDF; provided by Boost
    /// - Mode = 0; Median = 0; Mean = 0 for ν > 1 (undefined otherwise)
    /// - Variance = ν/(ν−2) for ν > 2 (∞ for 1 < ν ≤ 2; undefined for ν ≤ 1)
    ///
    /// This type constructs the underlying Boost `students_t_distribution` once
    /// and keeps an opaque handle; all evaluations reuse this instance.
    struct StudentT<T: BinaryFloatingPoint & Sendable>: Sendable, DistributionProtocol {
        typealias Real = T
        private final class Box: @unchecked Sendable {
            let raw: UnsafeRawPointer
            let free: (UnsafeMutableRawPointer) -> Void
            init(raw: UnsafeRawPointer, free: @escaping (UnsafeMutableRawPointer) -> Void) { self.raw = raw; self.free = free }
            deinit { free(UnsafeMutableRawPointer(mutating: raw)) }
        }

        public let degreesOfFreedom: T
        private let box: Box

        /// Initialize with ν degrees of freedom (> 0).
        public init(degreesOfFreedom v: T) throws {
            guard v > 0 else { throw DistributionError.parameterNotPositive(name: "degreesOfFreedom") }
            self.degreesOfFreedom = v
            if T.self == Double.self {
                guard let h = bs_student_t_make_d(Double(v)) else {
                    throw DistributionError.generalError(msg: "Distribution not initialized")
                }
                self.box = Box(raw: UnsafeRawPointer(h), free: { bs_student_t_free_d($0) })
            } else if T.self == Float.self {
                guard let h = bs_student_t_make_f(Float(v)) else {
                    throw DistributionError.generalError(msg: "Distribution not initialized")
                }
                self.box = Box(raw: UnsafeRawPointer(h), free: { bs_student_t_free_f($0) })
            } else {
                #if arch(x86_64) || arch(i386)
                guard let h = bs_student_t_make_l(Float80(v)) else {
                    throw DistributionError.generalError(msg: "Distribution not initialized")
                }
                self.box = Box(raw: UnsafeRawPointer(h), free: { bs_student_t_free_l($0) })
                #else
                guard let h = bs_student_t_make_d(Double(v)) else {
                    throw DistributionError.generalError(msg: "Distribution not initialized")
                }
                self.box = Box(raw: UnsafeRawPointer(h), free: { bs_student_t_free_d($0) })
                #endif
            }
        }

        // MARK: DistributionProtocol — Support
        public var supportLowerBound: T {
            if T.self == Double.self { let r = bs_student_t_range_d(box.raw); return T(r.lower) }
            if T.self == Float.self  { let r = bs_student_t_range_f(box.raw); return T(r.lower) }
            #if arch(x86_64) || arch(i386)
            let r = bs_student_t_range_l(box.raw); return T(r.lower)
            #else
            let r = bs_student_t_range_d(box.raw);  return T(r.lower)
            #endif
        }
        public var supportUpperBound: T {
            if T.self == Double.self { let r = bs_student_t_range_d(box.raw); return T(r.upper) }
            if T.self == Float.self  { let r = bs_student_t_range_f(box.raw); return T(r.upper) }
            #if arch(x86_64) || arch(i386)
            let r = bs_student_t_range_l(box.raw); return T(r.upper)
            #else
            let r = bs_student_t_range_d(box.raw);  return T(r.upper)
            #endif
        }
        public var range: (lower: T, upper: T) {
            if T.self == Double.self { let r = bs_student_t_range_d(box.raw);   return (T(r.lower), T(r.upper)) }
            if T.self == Float.self  { let r = bs_student_t_range_f(box.raw); return (T(r.lower), T(r.upper)) }
            #if arch(x86_64) || arch(i386)
            let r = bs_student_t_range_l(box.raw); return (T(r.lower), T(r.upper))
            #else
            let r = bs_student_t_range_d(box.raw);   return (T(r.lower), T(r.upper))
            #endif
        }

        // MARK: PDF/CDF/SF/Quantile
        public func pdf(_ x: T) throws -> T {
            if T.self == Double.self { return T(bs_student_t_pdf_d(box.raw, Double(x))) }
            if T.self == Float.self { return T(bs_student_t_pdf_f(box.raw, Float(x))) }
            #if arch(x86_64) || arch(i386)
            return T(bs_student_t_pdf_l(box.raw, Float80(x)))
            #else
            return T(bs_student_t_pdf_d(box.raw, Double(x)))
            #endif
        }

        public func logPdf(_ x: T) throws -> T {
            let p = try pdf(x)
            return T(log(Double(p)))
        }

        public func cdf(_ x: T) throws -> T {
            if T.self == Double.self { return T(bs_student_t_cdf_d(box.raw, Double(x))) }
            if T.self == Float.self { return T(bs_student_t_cdf_f(box.raw, Float(x))) }
            #if arch(x86_64) || arch(i386)
            return T(bs_student_t_cdf_l(box.raw, Float80(x)))
            #else
            return T(bs_student_t_cdf_d(box.raw, Double(x)))
            #endif
        }

        public func sf(_ x: T) throws -> T {
            if T.self == Double.self { return T(bs_student_t_ccdf_d(box.raw, Double(x))) }
            if T.self == Float.self { return T(bs_student_t_ccdf_f(box.raw, Float(x))) }
            #if arch(x86_64) || arch(i386)
            return T(bs_student_t_ccdf_l(box.raw, Float80(x)))
            #else
            return T(bs_student_t_ccdf_d(box.raw, Double(x)))
            #endif
        }

        public func quantile(_ p: T) throws -> T {
            if T.self == Double.self { return T(bs_student_t_quantile_d(box.raw, Double(p))) }
            if T.self == Float.self { return T(bs_student_t_quantile_f(box.raw, Float(p))) }
            #if arch(x86_64) || arch(i386)
            return T(bs_student_t_quantile_l(box.raw, Float80(p)))
            #else
            return T(bs_student_t_quantile_d(box.raw, Double(p)))
            #endif
        }

        public func quantileComplement(_ q: T) throws -> T {
            if T.self == Double.self { return T(bs_student_t_quantile_complement_d(box.raw, Double(q))) }
            if T.self == Float.self { return T(bs_student_t_quantile_complement_f(box.raw, Float(q))) }
            #if arch(x86_64) || arch(i386)
            return T(bs_student_t_quantile_complement_l(box.raw, Float80(q)))
            #else
            return T(bs_student_t_quantile_complement_d(box.raw, Double(q)))
            #endif
        }

        // MARK: Moments
        public var mean: T? { degreesOfFreedom > 1 ? T(0) : nil }
        public var variance: T? {
            let v: T = {
                if T.self == Double.self { return T(bs_student_t_variance_d(box.raw)) }
                if T.self == Float.self  { return T(bs_student_t_variance_f(box.raw)) }
                #if arch(x86_64) || arch(i386)
                return T(bs_student_t_variance_l(box.raw))
                #else
                return T(bs_student_t_variance_d(box.raw))
                #endif
            }()
            return v.isFinite ? v : nil
        }
        public var mode: T? { 0 }
        public var median: T { 0 }
        public var skewness: T? {
            let s: T = {
                if T.self == Double.self { return T(bs_student_t_skewness_d(box.raw)) }
                if T.self == Float.self  { return T(bs_student_t_skewness_f(box.raw)) }
                #if arch(x86_64) || arch(i386)
                return T(bs_student_t_skewness_l(box.raw))
                #else
                return T(bs_student_t_skewness_d(box.raw))
                #endif
            }()
            return s.isFinite ? s : nil
        }
        public var kurtosis: T? {
            let k: T = {
                if T.self == Double.self { return T(bs_student_t_kurtosis_d(box.raw)) }
                if T.self == Float.self  { return T(bs_student_t_kurtosis_f(box.raw)) }
                #if arch(x86_64) || arch(i386)
                return T(bs_student_t_kurtosis_l(box.raw))
                #else
                return T(bs_student_t_kurtosis_d(box.raw))
                #endif
            }()
            return k.isFinite ? k : nil
        }
        public var kurtosisExcess: T? {
            let e: T = {
                if T.self == Double.self { return T(bs_student_t_kurtosis_excess_d(box.raw)) }
                if T.self == Float.self  { return T(bs_student_t_kurtosis_excess_f(box.raw)) }
                #if arch(x86_64) || arch(i386)
                return T(bs_student_t_kurtosis_excess_l(box.raw))
                #else
                return T(bs_student_t_kurtosis_excess_d(box.raw))
                #endif
            }()
            return e.isFinite ? e : nil
        }

        // MARK: Hazards
        public func hazard(_ x: T) throws -> T {
            if T.self == Double.self { return T(bs_student_t_hazard_d(box.raw, Double(x))) }
            if T.self == Float.self  { return T(bs_student_t_hazard_f(box.raw, Float(x))) }
            #if arch(x86_64) || arch(i386)
            return T(bs_student_t_hazard_l(box.raw, Float80(x)))
            #else
            return T(bs_student_t_hazard_d(box.raw, Double(x)))
            #endif
        }
        public func chf(_ x: T) throws -> T {
            if T.self == Double.self { return T(bs_student_t_chf_d(box.raw, Double(x))) }
            if T.self == Float.self  { return T(bs_student_t_chf_f(box.raw, Float(x))) }
            #if arch(x86_64) || arch(i386)
            return T(bs_student_t_chf_l(box.raw, Float80(x)))
            #else
            return T(bs_student_t_chf_d(box.raw, Double(x)))
            #endif
        }
        
        // Lattice/discrete-only properties (continuous ⇒ nil)
        public var latticeStep: T? { nil }
        public var latticeOrigin: T? { nil }
        public var entropy: T? {
            if T.self == Double.self { return T(bs_student_t_entropy_d(box.raw)) }
            if T.self == Float.self  { return T(bs_student_t_entropy_f(box.raw)) }
            #if arch(x86_64) || arch(i386)
            return T(bs_student_t_entropy_l(box.raw))
            #else
            return T(bs_student_t_entropy_d(box.raw))
            #endif
        }


        // MARK: Planning helper
        /// Find degrees of freedom ν given effect size, α, β, and σ.
        public static func findDegreesOfFreedom(differenceFromMean: T, alpha: T, beta: T, sd: T, hint: T = 1) -> T {
            if T.self == Double.self {
                return T(bs_student_t_find_degrees_of_freedom_d(Double(differenceFromMean), Double(alpha), Double(beta), Double(sd), Double(hint)))
            }
            if T.self == Float.self {
                return T(bs_student_t_find_degrees_of_freedom_f(Float(differenceFromMean), Float(alpha), Float(beta), Float(sd), Float(hint)))
            }
            #if arch(x86_64) || arch(i386)
            return T(bs_student_t_find_degrees_of_freedom_l(Float80(differenceFromMean), Float80(alpha), Float80(beta), Float80(sd), Float80(hint)))
            #else
            return T(bs_student_t_find_degrees_of_freedom_d(Double(differenceFromMean), Double(alpha), Double(beta), Double(sd), Double(hint)))
            #endif
        }
    }
}
