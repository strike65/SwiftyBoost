//
//  Created by VT on 18.10.25.
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
    
    /// Gamma distribution Γ(k, θ) with shape k > 0 and scale θ > 0.
    ///
    /// - pdf(x) = x^(k−1) e^(−x/θ) / (Γ(k) θ^k), for x ≥ 0
    /// - cdf(x) = P(k, x/θ)  (regularized lower incomplete gamma)
    /// - sf(x)  = Q(k, x/θ)  (regularized upper incomplete gamma)
    /// - quantile(p) = θ * P⁻¹(k, p)
    ///
    /// All computations delegate to Boost via the SpecialFunctions wrappers.
    struct GammaDistribution<T: BinaryFloatingPoint> {
        // MARK: Identity
        public var description: String { "Gamma Distribution" }
        public var name: String { "gamma" }
        
        // MARK: Parameters
        public let shape: T  // k > 0
        public let scale: T  // θ > 0
        
        /// Initialize Γ(k, θ).
        ///
        /// - Parameters:
        ///   - shape: k > 0
        ///   - scale: θ > 0 (default 1)
        /// - Throws:
        ///   - SpecialFunctionError.parameterNotPositive if shape ≤ 0 or scale ≤ 0.
        public init(shape: T, scale: T = 1) throws {
            guard shape > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "shape") }
            guard scale > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "scale") }
            self.shape = shape
            self.scale = scale
        }
        
        // MARK: Moments and support
        /// Mean = k θ
        public var mean: T { shape * scale }
        /// Variance = k θ²
        public var variance: T { shape * scale * scale }
        /// Mode = (k − 1) θ for k ≥ 1; undefined for k < 1 (returns nil).
        public var mode: T? { shape >= 1 ? (shape - 1) * scale : nil }
        /// Support is x ∈ [0, +∞).
        public var supportLowerBound: T { 0 }
        public var supportUpperBound: T { T.infinity }
        
        // MARK: PDF / logPDF
        /// log pdf at x (throws if x < 0).
        ///
        /// log f(x) = (k−1) ln x − x/θ − ln Γ(k) − k ln θ
        public func logPdf(_ x: T) throws -> T {
            guard x >= 0 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: Double.infinity) }
            if x == 0 {
                // limit behavior: for k < 1, pdf diverges; for k == 1, pdf(0) = 1/θ; for k > 1, pdf(0) = 0
                // log-pdf at 0 is:
                // - +∞ when k < 1 (we return +∞)
                // - log(1/θ) when k == 1
                // - −∞ when k > 1
                if shape < 1 { return T(Double.infinity) }
                if shape == 1 { return -T(log(Double(scale))) }
                return -T(Double.infinity)
            }
            let k = Double(shape)
            let th = Double(scale)
            let xd = Double(x)
            // Use Boost-backed logGamma for stability:
//            let lgk = Double(try SpecialFunctions.logGamma(shape))
//            let val = (k - 1.0) * log(xd) - (xd / th) - lgk - k * log(th)
            let res = bs_gamma_pdf(Double(self.shape), Double.self(scale), Double(x))
            return T(res)
        }
        
        /// pdf at x (throws if x < 0).
        public func pdf(_ x: T) throws -> T {
            let lp = try logPdf(x)
            // exp on generic T: compute in Double and convert
            let v = exp(Double(lp))
            return T(v)
        }
        
        // MARK: CDF / SF
        /// CDF F(x) = P(k, x/θ), throws if x < 0.
        public func cdf(_ x: T) throws -> T {
            guard x >= 0 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: Double.infinity) }
            if x == 0 { return 0 }
//            return T(bs_gamma_cdf(Double(shape), Double(scale), Double(x)))
            return try SpecialFunctions.regularizedGammaP(shape, x: x / scale)
        }
        
        /// Survival function S(x) = Q(k, x/θ), throws if x < 0.
        public func sf(_ x: T) throws -> T {
            guard x >= 0 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: Double.infinity) }
            if x == 0 { return 1 }
            return try SpecialFunctions.regularizedGammaQ(shape, x: x / scale)
        }
        
        // MARK: Quantile / Inverse-CDF
        /// Quantile (inverse CDF). Returns x such that F(x) = p.
        ///
        /// - Parameter p: Probability in [0, 1].
        /// - Throws:
        ///   - SpecialFunctionError.parameterOutOfRange if p ∉ [0, 1].
        public func quantile(_ p: T) throws -> T {
            guard p >= 0 && p <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "p", min: 0.0, max: 1.0) }
            if p == 0 { return 0 }
            if p == 1 { return T.infinity }
            let xScaled = try SpecialFunctions.regularizedGammaPInv(shape, p: p)
            return scale * xScaled
        }
    }
}
