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
@testable import SwiftyBoost

@Suite("Cardinal B-splines basic behavior")
struct CardinalBSplinesTests {

    #if arch(x86_64)
    private func expectSpecialFunctionError(
        _ block: () throws -> Void
    ) {
        do {
            try block()
            #expect(Bool(false), "Expected SpecialFunctionError")
        } catch _ as SpecialFunctionError<Double> {
            // Expected type
        } catch _ as SpecialFunctionError<Float> {
            // Expected type
        } catch _ as SpecialFunctionError<Float80> {
            // Expected type (Float80 path)
        } catch let error {
            #expect(Bool(false), "Unexpected error: \(error)")
        }
    }
    #else
    private func expectSpecialFunctionError(
        _ block: () throws -> Void
    ) {
        do {
            try block()
            #expect(Bool(false), "Expected SpecialFunctionError")
        } catch _ as SpecialFunctionError<Double> {
            // Expected type
        } catch _ as SpecialFunctionError<Float> {
            // Expected type
        } catch let error {
            #expect(Bool(false), "Unexpected error: \(error)")
        }
    }
    #endif

    @Test("Value is finite and non-negative for x in [-1, 1] (Double/Float)")
    func inRangeEvaluations() throws {
        // Representative orders
        let orders = [1,2,3,4,5,6,7,8,9,11,12,13,14,15,16,17,18,19,20]
        let xs: [Double] = stride(from: -1.0, through: 1.0, by: 0.25).map { $0 }

        for n in orders {
            for x in xs {
                let vD: Double = try SpecialFunctions.cardinal_B_Spline(n, x)
                print("Oder: \(n), Spline: \(vD)")
                #expect(vD.isFinite)
                #expect(vD >= 0)

                let vF: Float = try SpecialFunctions.cardinal_B_Spline(n, Float(x))
                #expect(vF.isFinite)
                #expect(vF >= 0)
            }
        }
    }

    @Test("Symmetry: B is even, B' is odd, B'' is even (Double)")
    func symmetryDouble() throws {
        let orders = [3, 4, 6] // Derivatives require n ≥ 3
        let xs: [Double] = stride(from: 0.0, through: 1.0, by: 0.2).map { $0 }
        let tol = 1e-12

        for n in orders {
            for x in xs {
                let xp = x
                let xn = -x

                let bPos: Double = try SpecialFunctions.cardinal_B_Spline(n, xp)
                let bNeg: Double = try SpecialFunctions.cardinal_B_Spline(n, xn)
                #expect(abs(bPos - bNeg) <= tol)

                let bpPos: Double = try SpecialFunctions.cardinal_B_Spline_prime(n, xp)
                let bpNeg: Double = try SpecialFunctions.cardinal_B_Spline_prime(n, xn)
                #expect(abs(bpNeg + bpPos) <= tol)

                let bppPos: Double = try SpecialFunctions.cardinal_B_Spline_double_prime(n, xp)
                let bppNeg: Double = try SpecialFunctions.cardinal_B_Spline_double_prime(n, xn)
                #expect(abs(bppPos - bppNeg) <= tol)
            }

            // Additionally, odd derivative should be near zero at x = 0
            let bp0: Double = try SpecialFunctions.cardinal_B_Spline_prime(n, 0.0)
            #expect(abs(bp0) <= tol)
        }
    }

    @Test("Cross-precision consistency (Double vs Float)")
    func crossPrecisionConsistency() throws {
        let n = 5
        let xs: [Double] = stride(from: -1.0, through: 1.0, by: 0.1).map { $0 }

        for x in xs {
            let d: Double = try SpecialFunctions.cardinal_B_Spline(n, x)
            let f: Float = try SpecialFunctions.cardinal_B_Spline(n, Float(x))
            #expect(abs(d - Double(f)) <= 5e-6)

            let dp: Double = try SpecialFunctions.cardinal_B_Spline_prime(n, x)
            let fp: Float = try SpecialFunctions.cardinal_B_Spline_prime(n, Float(x))
            #expect(abs(dp - Double(fp)) <= 5e-5)

            let dpp: Double = try SpecialFunctions.cardinal_B_Spline_double_prime(n, x)
            let fpp: Float = try SpecialFunctions.cardinal_B_Spline_double_prime(n, Float(x))
            #expect(abs(dpp - Double(fpp)) <= 5e-4)
        }
    }

    @Test("x outside [-1, 1] throws SpecialFunctionError for all overloads")
    func outOfRangeXThrows() {
        let n = 4
        let badXsD: [Double] = [-1.000001, 1.000001]
        let badXsF: [Float] = [-1.0001, 1.0001]

        for x in badXsD {
            expectSpecialFunctionError {
                _ = try SpecialFunctions.cardinal_B_Spline(n, x)
            }
            expectSpecialFunctionError {
                _ = try SpecialFunctions.cardinal_B_Spline_prime(n, x)
            }
            expectSpecialFunctionError {
                _ = try SpecialFunctions.cardinal_B_Spline_double_prime(n, x)
            }
        }

        for x in badXsF {
            expectSpecialFunctionError {
                _ = try SpecialFunctions.cardinal_B_Spline(n, x)
            }
            expectSpecialFunctionError {
                _ = try SpecialFunctions.cardinal_B_Spline_prime(n, x)
            }
            expectSpecialFunctionError {
                _ = try SpecialFunctions.cardinal_B_Spline_double_prime(n, x)
            }
        }

        #if arch(i386) || arch(x86_64)
        let badXsL: [Float80] = [-1.0001, 1.0001]
        for x in badXsL {
            expectSpecialFunctionError {
                _ = try SpecialFunctions.cardinal_B_Spline(n, x)
            }
            expectSpecialFunctionError {
                _ = try SpecialFunctions.cardinal_B_Spline_prime(n, x)
            }
            expectSpecialFunctionError {
                _ = try SpecialFunctions.cardinal_B_Spline_double_prime(n, x)
            }
        }
        #endif
    }

    @Test("n < 3 throws for derivative functions")
    func nTooSmallThrows() {
        let badNs = [Int.min, -1, 0, 1, 2]
        let xD: Double = 0.0
        let xF: Float = 0.0

        for n in badNs {
            expectSpecialFunctionError {
                _ = try SpecialFunctions.cardinal_B_Spline_prime(n, xD)
            }
            expectSpecialFunctionError {
                _ = try SpecialFunctions.cardinal_B_Spline_double_prime(n, xD)
            }
            expectSpecialFunctionError {
                _ = try SpecialFunctions.cardinal_B_Spline_prime(n, xF)
            }
            expectSpecialFunctionError {
                _ = try SpecialFunctions.cardinal_B_Spline_double_prime(n, xF)
            }

            #if arch(i386) || arch(x86_64)
            let xL: Float80 = 0.0
            expectSpecialFunctionError {
                _ = try SpecialFunctions.cardinal_B_Spline_prime(n, xL)
            }
            expectSpecialFunctionError {
                _ = try SpecialFunctions.cardinal_B_Spline_double_prime(n, xL)
            }
            #endif
        }
    }
}
