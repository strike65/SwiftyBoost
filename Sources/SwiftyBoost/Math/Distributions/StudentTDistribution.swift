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
    public struct StudentT<T: Real & BinaryFloatingPoint & Sendable>: Sendable,
        DistributionProtocol
    {
        typealias RealType = T
        public let degreesOfFreedom: T
        private let dyn: Distribution.Dynamic<T>
        
        /// Initialize with ν degrees of freedom (> 0).
        public init(degreesOfFreedom v: T) throws {
            guard v > 0 else {
                throw DistributionError.parameterNotPositive(
                    name: "degreesOfFreedom"
                )
            }
            self.degreesOfFreedom = v
            self.dyn = try Distribution.Dynamic<T>(
                distributionName: "student_t",
                parameters: [
                    "df": v
                ]
            )
        }
        
        // MARK: DistributionProtocol — Support
        public var supportLowerBound: T { dyn.supportLowerBound }
        public var supportUpperBound: T { dyn.supportUpperBound }
        public var range: (lower: T, upper: T) { dyn.range }
        
        // MARK: PDF/CDF/SF/Quantile
        public func pdf(_ x: T) throws -> T { try dyn.pdf(x) }
        public func logPdf(_ x: T) throws -> T { try dyn.logPdf(x) }
        public func cdf(_ x: T) throws -> T { try dyn.cdf(x) }
        public func sf(_ x: T) throws -> T { try dyn.sf(x) }
        public func quantile(_ p: T) throws -> T { try dyn.quantile(p) }
        public func quantileComplement(_ q: T) throws -> T {
            try dyn.quantileComplement(q)
        }
        
        // MARK: Moments
        public var mean: T? { degreesOfFreedom > 1 ? T(0) : nil }
        public var variance: T? { dyn.variance }
        public var mode: T? { dyn.mode }
        public var median: T { dyn.median }
        public var skewness: T? { dyn.skewness }
        public var kurtosis: T? { dyn.kurtosis }
        public var kurtosisExcess: T? { dyn.kurtosisExcess }
        
        // MARK: Hazards
        public func hazard(_ x: T) throws -> T { try dyn.hazard(x) }
        public func chf(_ x: T) throws -> T { try dyn.chf(x) }
        
        // Lattice/discrete-only properties (continuous ⇒ nil)
        public var latticeStep: T? { nil }
        public var latticeOrigin: T? { nil }
        public var entropy: T? { dyn.entropy }
        
        // MARK: Planning helper
        /// Find degrees of freedom ν given effect size, α, β, and σ.
        public static func findDegreesOfFreedom(
            differenceFromMean: T,
            alpha: T,
            beta: T,
            sd: T,
            hint: T = 1
        ) -> Int {
            let res: T
            if T.self == Double.self {
                res = T(
                    ceil(bs_student_t_find_degrees_of_freedom_d(
                        Double(differenceFromMean),
                        Double(alpha),
                        Double(beta),
                        Double(sd),
                        Double(hint)
                    ))
                )
            }
            else if T.self == Float.self {
                res = T(
                    ceil(bs_student_t_find_degrees_of_freedom_f(
                        Float(differenceFromMean),
                        Float(alpha),
                        Float(beta),
                        Float(sd),
                        Float(hint)
                    )))
            }
            else {
#if arch(x86_64) || arch(i386)
                res = T(
                    ceil(bs_student_t_find_degrees_of_freedom_l(
                        Float80(differenceFromMean),
                        Float80(alpha),
                        Float80(beta),
                        Float80(sd),
                        Float80(hint)
                    )))
#else
                res = T(
                    ceil(bs_student_t_find_degrees_of_freedom_d(
                        Double(differenceFromMean),
                        Double(alpha),
                        Double(beta),
                        Double(sd),
                        Double(hint)
                    )))
#endif
            }
            return Int(res) + 1
        }
    }
}
