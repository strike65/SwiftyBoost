//
//  Created by VT on 17.10.25.
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

import Testing
import Darwin
@testable import SwiftyBoost

@Suite("Gegenbauer polynomial tests (Double)")
struct GegenbauerDoubleTests {

    private func expectSpecialFunctionError(
        _ block: () throws -> Void
    ) {
        do {
            try block()
            #expect(Bool(false), "Expected SpecialFunctionError")
        } catch _ as SpecialFunctionError<Double> {
            // Expected
        } catch {
            #expect(Bool(false), "Unexpected error: \(error)")
        }
    }

    // Tolerance for Double
    let tol = 1e-12

    // Small helpers for known polynomials
    // Chebyshev U_n(x) for n = 0...3
    func chebyshevU(_ n: Int, _ x: Double) -> Double {
        switch n {
        case 0: return 1
        case 1: return 2*x
        case 2: return 4*x*x - 1
        case 3: return 8*x*x*x - 4*x
        default:
            // 3-term recurrence for U_n: U_n = 2x U_{n-1} - U_{n-2}
            var u0 = 1.0
            var u1 = 2.0 * x
            if n == 0 { return u0 }
            if n == 1 { return u1 }
            for _ in 2...n {
                let u = 2.0 * x * u1 - u0
                u0 = u1
                u1 = u
            }
            return u1
        }
    }

    // Legendre P_n(x) for n = 0...3
    func legendreP(_ n: Int, _ x: Double) -> Double {
        switch n {
        case 0: return 1
        case 1: return x
        case 2: return 0.5 * (3*x*x - 1)
        case 3: return 0.5 * (5*x*x*x - 3*x)
        default:
            // standard recurrence: (n+1) P_{n+1} = (2n+1) x P_n - n P_{n-1}
            var p0 = 1.0
            var p1 = x
            if n == 0 { return p0 }
            if n == 1 { return p1 }
            for k in 1..<(n) {
                let pkp1 = ((2.0*Double(k)+1.0) * x * p1 - Double(k) * p0) / Double(k+1)
                p0 = p1
                p1 = pkp1
            }
            return p1
        }
    }

    // Rising Pochhammer (λ)_k
    func pochhammer(_ a: Double, _ k: Int) -> Double {
        if k <= 0 { return 1 }
        var prod = 1.0
        for i in 0..<k {
            prod *= (a + Double(i))
        }
        return prod
    }

    @Test("Base values: C0^λ(x) = 1, C1^λ(x) = 2λx, various λ and x")
    func baseValues() throws {
        let lambdas: [Double] = [0.0, 0.5, 1.0, 2.0]
        let xs: [Double] = [-0.9, -0.2, 0.0, 0.3, 0.8]
        for λ in lambdas {
            for x in xs {
                let c0 = try SpecialFunctions.gegenbauer(n: 0, lambda: λ, x: x)
                #expect(abs(c0 - 1.0) <= tol, "C0 should be 1")

                let c1 = try SpecialFunctions.gegenbauer(n: 1, lambda: λ, x: x)
                let expected = 2.0 * λ * x
                #expect(abs(c1 - expected) <= 1e-12, "C1 should be 2λx")
            }
        }
    }

    @Test("Special case λ = 0: C0^0 = 1 and Cn^0 = 0 for n > 0")
    func lambdaZero() throws {
        let xs: [Double] = [-0.9, -0.2, 0.0, 0.3, 0.8]
        for x in xs {
            let c0 = try SpecialFunctions.gegenbauer(n: 0, lambda: 0.0, x: x)
            #expect(abs(c0 - 1.0) <= tol)

            for n in 1...6 {
                let cn = try SpecialFunctions.gegenbauer(n: n, lambda: 0.0, x: x)
                #expect(abs(cn) <= tol, "Cn^0 should be 0 for n>0")
            }
        }
    }

    @Test("Chebyshev identity: Cn^1(x) = U_n(x) for small n")
    func chebyshevUIdentity() throws {
        let xs: [Double] = [-0.9, -0.2, 0.0, 0.3, 0.8]
        for n in 0...8 {
            for x in xs {
                let cn = try SpecialFunctions.gegenbauer(n: n, lambda: 1.0, x: x)
                let un = chebyshevU(n, x)
                #expect(abs(cn - un) <= 1e-12, "Mismatch for n=\(n), x=\(x)")
            }
        }
    }

    @Test("Legendre identity: Cn^(1/2)(x) = P_n(x) for small n")
    func legendreIdentity() throws {
        let xs: [Double] = [-0.9, -0.2, 0.0, 0.3, 0.8]
        for n in 0...6 {
            for x in xs {
                let cn = try SpecialFunctions.gegenbauer(n: n, lambda: 0.5, x: x)
                let pn = legendreP(n, x)
                #expect(abs(cn - pn) <= 5e-13, "Mismatch for n=\(n), x=\(x)")
            }
        }
    }

    @Test("First derivative identity: d/dx Cn^λ(x) = 2λ C_{n-1}^{λ+1}(x)")
    func firstDerivativeIdentity() throws {
        let xs: [Double] = [-0.8, -0.1, 0.0, 0.4, 0.9]
        let lambdas: [Double] = [0.5, 1.0, 2.5]
        for λ in lambdas {
            for n in 1...8 {
                for x in xs {
                    let lhs = try SpecialFunctions.gegenbauerPrime(n: n, lambda: λ, x: x)
                    let rhsPoly = try SpecialFunctions.gegenbauer(n: n-1, lambda: λ + 1.0, x: x)
                    let rhs = 2.0 * λ * rhsPoly
                    #expect(abs(lhs - rhs) <= 2e-12, "n=\(n), λ=\(λ), x=\(x)")
                }
            }
        }
    }

    @Test("k-th derivative identity: d^k/dx^k Cn^λ = 2^k (λ)_k C_{n-k}^{λ+k}, for k ≤ n")
    func kthDerivativeIdentity() throws {
        let xs: [Double] = [-0.7, -0.2, 0.0, 0.5]
        let lambdas: [Double] = [0.5, 1.0, 1.75]
        for λ in lambdas {
            for n in 2...8 {
                for k in 0...min(3, n) {
                    for x in xs {
                        let lhs = try SpecialFunctions.gegenbauerDerivative(n: n, lambda: λ, x: x, k: k)
                        let coeff = Darwin.pow(2.0, Double(k)) * pochhammer(λ, k)
                        let poly = try SpecialFunctions.gegenbauer(n: n - k, lambda: λ + Double(k), x: x)
                        let rhs = coeff * poly
                        #expect(abs(lhs - rhs) <= 5e-12, "n=\(n), k=\(k), λ=\(λ), x=\(x)")
                    }
                }
            }
        }
    }

    @Test("Error: n < 0 throws parameterNotPositive")
    func errorNegativeN() async {
        expectSpecialFunctionError {
            _ = try SpecialFunctions.gegenbauer(n: -1, lambda: 1.0, x: 0.2)
        }
    }

    @Test("Error: k < 0 throws parameterNotPositive")
    func errorNegativeK() async {
        expectSpecialFunctionError {
            _ = try SpecialFunctions.gegenbauerDerivative(n: 3, lambda: 1.0, x: 0.2, k: -1)
        }
    }

    @Test("Error: non-finite lambda and x throw parameterNotFinite")
    func errorNonFinite() async {
        expectSpecialFunctionError {
            _ = try SpecialFunctions.gegenbauer(n: 2, lambda: Double.nan, x: 0.0)
        }

        expectSpecialFunctionError {
            _ = try SpecialFunctions.gegenbauer(n: 2, lambda: 1.0, x: Double.infinity)
        }
    }
}

@Suite("Gegenbauer polynomial tests (Float)")
struct GegenbauerFloatTests {
    let tol: Float = 1e-5

    @Test("Chebyshev identity for Float: Cn^1(x) = U_n(x)")
    func chebyshevFloat() throws {
        func U(_ n: Int, _ x: Float) -> Float {
            switch n {
            case 0: return 1
            case 1: return 2*x
            case 2: return 4*x*x - 1
            case 3: return 8*x*x*x - 4*x
            default:
                var u0: Float = 1
                var u1: Float = 2*x
                if n == 0 { return u0 }
                if n == 1 { return u1 }
                for _ in 2...n {
                    let u = 2*x*u1 - u0
                    u0 = u1
                    u1 = u
                }
                return u1
            }
        }

        let xs: [Float] = [-0.8, -0.1, 0.0, 0.4, 0.9]
        for n in 0...8 {
            for x in xs {
                let cn = try SpecialFunctions.gegenbauer(n: n, lambda: Float(1), x: x)
                let un = U(n, x)
                #expect(abs(cn - un) <= tol, "n=\(n), x=\(x)")
            }
        }
    }

    @Test("First derivative identity for Float")
    func firstDerivativeFloat() throws {
        let xs: [Float] = [-0.7, 0.0, 0.6]
        for n in 1...6 {
            for x in xs {
                let λ: Float = 1.25
                let lhs = try SpecialFunctions.gegenbauerPrime(n: n, lambda: λ, x: x)
                let rhsPoly = try SpecialFunctions.gegenbauer(n: n-1, lambda: λ + 1, x: x)
                let rhs = 2*λ*rhsPoly
                #expect(abs(lhs - rhs) <= 5e-5, "n=\(n), x=\(x)")
            }
        }
    }
}

#if arch(x86_64)
@Suite("Gegenbauer polynomial tests (Float80)")
struct GegenbauerFloat80Tests {
    let tol: Float80 = 1e-15

    @Test("Legendre identity for Float80: Cn^(1/2)(x) = P_n(x)")
    func legendreFloat80() throws {
        func P(_ n: Int, _ x: Float80) -> Float80 {
            switch n {
            case 0: return 1
            case 1: return x
            case 2: return 0.5 * (3*x*x - 1)
            case 3: return 0.5 * (5*x*x*x - 3*x)
            default:
                var p0: Float80 = 1
                var p1: Float80 = x
                if n == 0 { return p0 }
                if n == 1 { return p1 }
                for k in 1..<(n) {
                    let pkp1 = ((2*Float80(k)+1) * x * p1 - Float80(k) * p0) / Float80(k+1)
                    p0 = p1
                    p1 = pkp1
                }
                return p1
            }
        }

        let xs: [Float80] = [-0.8, -0.1, 0.0, 0.4, 0.9]
        for n in 0...8 {
            for x in xs {
                let cn = try SpecialFunctions.gegenbauer(n: n, lambda: Float80(0.5), x: x)
                let pn = P(n, x)
                #expect(abs(cn - pn) <= tol, "n=\(n), x=\(x)")
            }
        }
    }
}
#endif
