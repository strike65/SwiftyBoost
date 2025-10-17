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

@Suite("Bessel J zeros")
struct BesselZerosTests {

    // Helper to assert approximate equality
    func expectApproximatelyEqual(_ a: Double, _ b: Double, tol: Double, _ message: String = "") {
        #expect(abs(a - b) <= tol, "\(message) | got \(a), expected \(b), tol \(tol)")
    }

    @Test("J0 first zero")
    func j0FirstZero() async throws {
        // Reference: J0 first zero ≈ 2.404825557695772768621...
        let z: Double = SpecialFunctions.besselJZero(v: 0.0, m: 1)
        expectApproximatelyEqual(z, 2.4048255576957728, tol: 1e-13, "J0 first zero")
    }

    @Test("J2 first five zeros")
    func j2FirstFiveZeros() async throws {
        // Reference values (MathWorld / Boost example)
        let ref: [Double] = [
            5.13562230184068,
            8.41724414039986,
            11.6198411721491,
            14.7959517823513,
            17.9598194949878
        ]
        let zs = SpecialFunctions.besselJZeros(v: 2.0, start: 1, count: ref.count)
        #expect(zs.count == ref.count)

        for (i, (z, r)) in zip(zs, ref).enumerated() {
            expectApproximatelyEqual(Double(z), r, tol: 5e-13, "J2 zero #\(i + 1)")
        }
    }

    @Test("Fractional order v = 71/19, first three zeros")
    func fractionalOrderFirstThree() async throws {
        // Reference values from Boost example / Mathematica (50 digits trimmed)
        let v = 71.0 / 19.0
        let ref: [Double] = [
            7.2731751938316489503185694262290765588963196701623,
            10.724858308883141732536172745851416647110749599085,
            14.018504599452388106120459558042660282427471931581
        ]
        let zs = SpecialFunctions.besselJZeros(v: v, start: 1, count: 3)
        #expect(zs.count == 3)

        // Double precision can’t match 50 digits; 1e-12 is a reasonable tolerance here.
        for (i, (z, r)) in zip(zs, ref).enumerated() {
            expectApproximatelyEqual(Double(z), r, tol: 1e-12, "J_{\(v)} zero #\(i + 1)")
        }
    }

    @Test("Empty result for count = 0")
    func emptyCount() async throws {
        let zs = SpecialFunctions.besselJZeros(v: 0.0, start: 1, count: 0)
        #expect(zs.isEmpty)
    }

    @Test("Start index offset behavior")
    func startIndexOffset() async throws {
        // Compare first three zeros to the next three when starting at 2
        let firstThree = SpecialFunctions.besselJZeros(v: 0.0, start: 1, count: 3)
        let nextTwo = SpecialFunctions.besselJZeros(v: 0.0, start: 2, count: 2)

        #expect(firstThree.count == 3)
        #expect(nextTwo.count == 2)

        // The last two of the firstThree should match the two from start:2
        expectApproximatelyEqual(Double(firstThree[1]), Double(nextTwo[0]), tol: 1e-13, "Start index consistency (2nd vs 1st)")
        expectApproximatelyEqual(Double(firstThree[2]), Double(nextTwo[1]), tol: 1e-13, "Start index consistency (3rd vs 2nd)")
    }

    // MARK: - Float coverage

    @Test("J0 first zero (Float)")
    func j0FirstZeroFloat() async throws {
        let z: Float = SpecialFunctions.besselJZero(v: 0 as Float, m: 1)
        expectApproximatelyEqual(Double(z), 2.4048255576957728, tol: 2e-6, "J0 first zero (Float)")
    }

    @Test("J2 first five zeros (Float)")
    func j2FirstFiveZerosFloat() async throws {
        let ref: [Double] = [
            5.13562230184068,
            8.41724414039986,
            11.6198411721491,
            14.7959517823513,
            17.9598194949878
        ]
        let zs: [Float] = SpecialFunctions.besselJZeros(v: 2 as Float, start: 1, count: ref.count)
        #expect(zs.count == ref.count)
        for (i, (z, r)) in zip(zs, ref).enumerated() {
            expectApproximatelyEqual(Double(z), r, tol: 2e-6, "J2 zero (Float) #\(i + 1)")
        }
    }

    @Test("Fractional order v = 71/19, first three zeros (Float)")
    func fractionalOrderFirstThreeFloat() async throws {
        let vF = Float(71.0 / 19.0)
        let ref: [Double] = [
            7.2731751938316489503185694262290765588963196701623,
            10.724858308883141732536172745851416647110749599085,
            14.018504599452388106120459558042660282427471931581
        ]
        let zs: [Float] = SpecialFunctions.besselJZeros(v: vF, start: 1, count: 3)
        #expect(zs.count == 3)
        for (i, (z, r)) in zip(zs, ref).enumerated() {
            expectApproximatelyEqual(Double(z), r, tol: 2e-6, "J_{\(Double(vF))} zero (Float) #\(i + 1)")
        }
    }

    @Test("Empty result for count = 0 (Float)")
    func emptyCountFloat() async throws {
        let zs: [Float] = SpecialFunctions.besselJZeros(v: 0 as Float, start: 1, count: 0)
        #expect(zs.isEmpty)
    }

    @Test("Start index offset behavior (Float)")
    func startIndexOffsetFloat() async throws {
        let firstThree: [Float] = SpecialFunctions.besselJZeros(v: 0 as Float, start: 1, count: 3)
        let nextTwo: [Float] = SpecialFunctions.besselJZeros(v: 0 as Float, start: 2, count: 2)
        #expect(firstThree.count == 3)
        #expect(nextTwo.count == 2)
        expectApproximatelyEqual(Double(firstThree[1]), Double(nextTwo[0]), tol: 2e-6, "Start index consistency (Float) (2nd vs 1st)")
        expectApproximatelyEqual(Double(firstThree[2]), Double(nextTwo[1]), tol: 2e-6, "Start index consistency (Float) (3rd vs 2nd)")
    }

    // MARK: - Float80 coverage (x86_64)

    #if arch(x86_64)
    @Test("J0 first zero (Float80)")
    func j0FirstZeroFloat80() async throws {
        let z: Float80 = SpecialFunctions.besselJZero(v: 0 as Float80, m: 1)
        expectApproximatelyEqual(Double(z), 2.4048255576957728, tol: 1e-18, "J0 first zero (Float80)")
    }

    @Test("J2 first five zeros (Float80)")
    func j2FirstFiveZerosFloat80() async throws {
        let ref: [Double] = [
            5.13562230184068,
            8.41724414039986,
            11.6198411721491,
            14.7959517823513,
            17.9598194949878
        ]
        let zs: [Float80] = SpecialFunctions.besselJZeros(v: 2 as Float80, start: 1, count: ref.count)
        #expect(zs.count == ref.count)
        for (i, (z, r)) in zip(zs, ref).enumerated() {
            expectApproximatelyEqual(Double(z), r, tol: 1e-18, "J2 zero (Float80) #\(i + 1)")
        }
    }

    @Test("Fractional order v = 71/19, first three zeros (Float80)")
    func fractionalOrderFirstThreeFloat80() async throws {
        let vL = Float80(71.0) / 19.0
        let ref: [Double] = [
            7.2731751938316489503185694262290765588963196701623,
            10.724858308883141732536172745851416647110749599085,
            14.018504599452388106120459558042660282427471931581
        ]
        let zs: [Float80] = SpecialFunctions.besselJZeros(v: vL, start: 1, count: 3)
        #expect(zs.count == 3)
        for (i, (z, r)) in zip(zs, ref).enumerated() {
            expectApproximatelyEqual(Double(z), r, tol: 1e-18, "J_{\(Double(vL))} zero (Float80) #\(i + 1)")
        }
    }

    @Test("Empty result for count = 0 (Float80)")
    func emptyCountFloat80() async throws {
        let zs: [Float80] = SpecialFunctions.besselJZeros(v: 0 as Float80, start: 1, count: 0)
        #expect(zs.isEmpty)
    }

    @Test("Start index offset behavior (Float80)")
    func startIndexOffsetFloat80() async throws {
        let firstThree: [Float80] = SpecialFunctions.besselJZeros(v: 0 as Float80, start: 1, count: 3)
        let nextTwo: [Float80] = SpecialFunctions.besselJZeros(v: 0 as Float80, start: 2, count: 2)
        #expect(firstThree.count == 3)
        #expect(nextTwo.count == 2)
        expectApproximatelyEqual(Double(firstThree[1]), Double(nextTwo[0]), tol: 1e-18, "Start index consistency (Float80) (2nd vs 1st)")
        expectApproximatelyEqual(Double(firstThree[2]), Double(nextTwo[1]), tol: 1e-18, "Start index consistency (Float80) (3rd vs 2nd)")
    }
    #endif

}
