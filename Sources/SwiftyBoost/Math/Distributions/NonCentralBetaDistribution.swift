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
    public struct NonCentralBeta<T: Real & BinaryFloatingPoint & Sendable>: Sendable, DistributionProtocol {
        typealias RealType = T
        public let alpha: T
        public let beta: T
        public let lambda: T
        private let dyn: Distribution.Dynamic<T>
        public init(alpha: T = 0, beta: T = 1, lambda: T) throws {
            guard beta > 0 else {
                throw DistributionError.parameterNotPositive(
                    name: "beta",
                    value: beta
                )
            }
            guard alpha > 0 else {
                throw DistributionError.parameterNotPositive(
                    name: "alpha",
                    value: alpha
                )
            }
            guard lambda >= 0 else {
                throw DistributionError.parameterOutOfRange(name: "lambda", min: 0, max: .infinity)
            }
            self.alpha = alpha
            self.beta = beta
            self.lambda = lambda
            self.dyn = try Distribution.Dynamic<T>(
                distributionName: "non_central_beta",
                parameters: [
                    "alpha": alpha,
                    "beta": beta,
                    "lambda": lambda
                ]
            )
        }
        public var supportLowerBound: T { dyn.supportLowerBound }
        public var supportUpperBound: T { dyn.supportUpperBound }
        public var range: (lower: T, upper: T) { dyn.range }
        public func pdf(_ x: T) throws -> T { try dyn.pdf(x) }
        public func logPdf(_ x: T) throws -> T { try dyn.logPdf(x) }
        public func cdf(_ x: T) throws -> T { try dyn.cdf(x) }
        public func sf(_ x: T) throws -> T { try dyn.sf(x) }
        public func quantile(_ p: T) throws -> T { try dyn.quantile(p) }
        public func quantileComplement(_ q: T) throws -> T {
            try dyn.quantileComplement(q)
        }
        public var mean: T? { nil }
        public var variance: T? { dyn.variance }
        public var mode: T? { dyn.mode }
        public var median: T { dyn.median }
        public var skewness: T? {
            guard self.moments.mu2.isFinite, self.moments.mu3.isFinite, self.moments.mu2 > 0 else { return nil }
            return self.moments.mu3 / T.pow(self.moments.mu2, T(1.5))
        }
        public var kurtosis: T? { self.moments.mu4 / (self.moments.mu2 * self.moments.mu2) }
        public var kurtosisExcess: T? { self.moments.mu4 / (self.moments.mu2 * self.moments.mu2) - 3 }
        public func hazard(_ x: T) throws -> T { try dyn.hazard(x) }
        public func chf(_ x: T) throws -> T { try dyn.chf(x) }
        public var latticeStep: T? { nil }
        public var latticeOrigin: T? { nil }
        public var entropy: T? {
            // split up complex expression (Swift type checker)
            guard self.alpha > 0, self.beta > 0, self.lambda >= 0 else { return nil }
            let one: T = 1
            let two: T = 2
            if self.lambda == 0 {
                do {
                    let logB1 = try SpecialFunctions.logGamma(self.alpha)
                    let logB2 = try SpecialFunctions.logGamma(self.beta)
                    let logAb2 = try SpecialFunctions.logGamma(self.alpha + self.beta)
                    let logB = logB1 + logB2 - logAb2

                    let dA = try SpecialFunctions.digamma(self.alpha)
                    let dB = try SpecialFunctions.digamma(self.beta)
                    let dAB = try SpecialFunctions.digamma(self.alpha + self.beta)

                    let term1 = -(self.alpha - one) * dA
                    let term2 = term1 - (self.beta - one) * dB
                    let term3 = term2 + (self.alpha + self.beta - two) * dAB

                    let h = logB + term3
                    return h
                } catch { return nil }
            }
            do {
                let pTail: T = 1e-15
                let tol:   T = 1e-10
                let mu = self.lambda / two
                guard mu > 0 else { return nil }
                let (Jm, Jp) = try poisson_window(mean: mu, tail: pTail)
                let logmu: T = T.log(mu)
                let lgB: T   = try SpecialFunctions.logGamma(self.beta)
                let logw: [T] = try (Jm...Jp).map { j in
                    let jj  = T(j)
                    let lg1: T = try SpecialFunctions.logGamma(jj + one)
                    return -mu + jj * logmu - lg1
                }
                let lnbet: [T] = try (Jm...Jp).map { j in
                    let jj   = T(j)
                    let lgA: T  = try SpecialFunctions.logGamma(self.alpha + jj)
                    let lgAB: T = try SpecialFunctions.logGamma(self.alpha + self.beta + jj)
                    return lgA + lgB - lgAB
                }
                let integrand: @Sendable (T) -> T = { (x: T) in
                    let xx = min(max(x, T.leastNonzeroMagnitude), one - T.ulpOfOne)
                    let lx: T = T.log(xx)
                    let l1x: T = T.log(onePlus: -xx)
                    
                    var m: T = -T.infinity
                    for j in Jm...Jp {
                        let idx: Int  = j - Jm
                        let jj: T = T(j)
                        let term1: T = (self.alpha + jj - one) * lx
                        let term2: T = (self.beta  - one) * l1x
                        let term: T = logw[idx] + term1 + term2 - lnbet[idx]
                        if term > m { m = term }
                    }
                    var s: T = 0
                    for j in Jm...Jp {
                        let idx: Int = j - Jm
                        let jj: T = T(j)
                        let term1: T = (self.alpha + jj - one) * lx
                        let term2: T = (self.beta  - one) * l1x
                        let term = logw[idx] + term1 + term2 - lnbet[idx]
                        s += T.exp(term - m)
                    }
                    let logf: T = m + T.log(s)
                    let f: T = T.exp(logf)
                    return -f * logf
                }
                
                let integrator = try Quadrature.Integrator<T>(
                    rule: .tanhSinh(maxRefinements: 1000, tolerance: Double(tol))
                )
                let I: T = try integrator
                    .integrate(over: .finite(lower: 0, upper: 1), integrand: integrand)
                    .value
                
                return I
            } catch {
                return nil
            }
        }
        public var isDiscrete: Bool { false }
        public func klDivergence(
            relativeTo other: Self,
            options: Distribution.KLDivergenceOptions<T> = .automatic()
        ) throws -> T? {
            try dyn.klDivergence(relativeTo: other.dyn, options: options)
        }
        // moments
        public var moments: (mu1: T, mu2: T, mu3: T, mu4: T) {
            do {
                let one: T = 1
                let two: T = 2
                let three: T = 3
                let four: T = 4
                let six: T = 6
                let a:T = self.alpha
                let b:T = self.beta
                
                @inline(__always)
                func ratioRising(_ base: T, _ denBase: T, _ r: Int) -> T {
                    var num: T = 1, den: T = 1
                    if r >= 1 { num *= base; den *= denBase }
                    if r >= 2 { num *= (base + 1); den *= (denBase + 1) }
                    if r >= 3 { num *= (base + 2); den *= (denBase + 2) }
                    if r >= 4 { num *= (base + 3); den *= (denBase + 3) }
                    return num / den
                }
                
                if self.lambda.isZero {
                    let M11 = try SpecialFunctions.pochhammer(a, 1)
                    let M12 = try SpecialFunctions.pochhammer(a + b, 1)
                    let M1 = M11 / M12
                    let M21 = try SpecialFunctions.pochhammer(a, 2)
                    let M22 = try SpecialFunctions.pochhammer(a + b, 2)
                    let M2 = M21 / M22
                    let M31 = try SpecialFunctions.pochhammer(a, 3)
                    let M32 = try SpecialFunctions.pochhammer(a + b, 3)
                    let M3 = M31 / M32
                    let M41 = try SpecialFunctions.pochhammer(a, 4)
                    let M42 = try SpecialFunctions.pochhammer(a + b, 4)
                    let M4 = M41 / M42
                    let mu1 = M1
                    let mu2 = M2 - mu1 * mu1
                    let mu3 = M3 - three * mu1 * M2 + two * mu1 * mu1 * mu1
                    var mu4 = M4 - four * mu1 * M3
                    mu4 = mu4 + six * mu1 * mu1 * M2
                    mu4 = mu4 - three * mu1 * mu1 * mu1 * mu1
                    return (mu1, mu2, mu3, mu4)
                }
                let mu = self.lambda / two
                let (Jm,Jp) = try poisson_window(mean: mu, tail: T(1e-18))
                let logmu: T = T.log(mu)
                let logw: [T] = try (Jm...Jp).map { j in
                    let jj  = T(j)
                    let lg1: T = try SpecialFunctions.logGamma(jj + one)
                    return -mu + jj * logmu - lg1
                }
                var Mr: T = 0, c:T = 0
                var M:[T] = [T](repeating: T.zero, count: 5)
                let w = logw.map(T.exp)
                for r in 1...4 {
                    Mr = 0
                    c = 0
                    for j in Jm...Jp {
                        let idx = j - Jm
                        let jj = T(j)
                        let frac:T = ratioRising(a + jj, a + b + jj, r)
                        let y = w[idx] * frac - c
                        let t = Mr + y
                        c = (t - Mr) - y
                        Mr = t
                    }
                    M[r] = Mr
                }
                let mu1 = M[1]
                let mu2 = M[2] - mu1 * mu1
                let mu3 = M[3] - three * mu1 * M[2] + two * mu1 * mu1 * mu1
                var mu4 = M[4] - four * mu1 * M[3]
                mu4 = mu4 + six * mu1 * mu1 * M[2]
                mu4 = mu4 - three * mu1 * mu1 * mu1 * mu1
                return (mu1, mu2, mu3, mu4)
            }
            catch _ {
                return (.nan, .nan,.nan,.nan)
            }
        }
    }
}

