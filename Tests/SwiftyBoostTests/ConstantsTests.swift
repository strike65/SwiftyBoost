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

// Generic helpers
@inline(__always)
func approxEqual<T: BinaryFloatingPoint>(_ a: T, _ b: T, tol: T) -> Bool {
    let diff = a > b ? a - b : b - a
    return diff <= tol
}

// MARK: - Double tests

@Suite("Constants<Double> tests")
struct DoubleConstantsTests {
    let tol: Double = 1e-15

    @Test
    func piFamilyAndDerived() {
        let pi = Constants<Double>.pi
        #expect(pi.isFinite && pi > 3 && pi < 4)

        #expect(approxEqual(Constants<Double>.twoPi, 2 * pi, tol: tol))
        #expect(approxEqual(Constants<Double>.halfPi, pi / 2, tol: tol))
        #expect(approxEqual(Constants<Double>.quarterPi, pi / 4, tol: tol))
        #expect(approxEqual(Constants<Double>.thirdPi, pi / 3, tol: tol))
        #expect(approxEqual(Constants<Double>.twoThirdsPi, 2 * pi / 3, tol: tol))
        #expect(approxEqual(Constants<Double>.threeQuartersPi, 3 * pi / 4, tol: tol))
        #expect(approxEqual(Constants<Double>.sixthPi, pi / 6, tol: tol))
        #expect(approxEqual(Constants<Double>.piSqr, pi * pi, tol: tol))
    }

    @Test
    func rootsAndReciprocals() {
        let pi = Constants<Double>.pi
        let rootTwo = Constants<Double>.rootTwo
        let rootThree = Constants<Double>.rootThree
        let rootPi = Constants<Double>.rootPi

        #expect(approxEqual(rootTwo * rootTwo, 2.0, tol: tol))
        #expect(approxEqual(rootThree * rootThree, 3.0, tol: tol))
        #expect(approxEqual(rootPi * rootPi, pi, tol: tol))

        #expect(approxEqual(Constants<Double>.oneDivPi, 1.0 / pi, tol: tol))
        #expect(approxEqual(Constants<Double>.oneDivTwoPi, 1.0 / (2.0 * pi), tol: tol))
        #expect(approxEqual(Constants<Double>.oneDivRootPi, 1.0 / rootPi, tol: tol))
        #expect(approxEqual(Constants<Double>.twoDivPi, 2.0 / pi, tol: tol))
        #expect(approxEqual(Constants<Double>.twoDivRootPi, 2.0 / rootPi, tol: tol))
    }

    @Test
    func logsAndExponential() {
        let ln2 = Constants<Double>.lnTwo
        let ln10 = Constants<Double>.lnTen
        let lnln2 = Constants<Double>.lnLnTwo
        let e = Constants<Double>.e

        #expect(approxEqual(ln2, Darwin.log(2.0), tol: tol))
        #expect(approxEqual(ln10, Darwin.log(10.0), tol: tol))
        #expect(approxEqual(lnln2, Darwin.log(ln2), tol: tol))
        #expect(approxEqual(e, Darwin.exp(1.0), tol: tol))
    }

    @Test
    func specialConstants() {
        // Reference values
        let gammaRef = 0.5772156649015328606
        let catalanRef = 0.91596559417721901505
        let zeta3Ref = 1.2020569031595942854

        #expect(approxEqual(Constants<Double>.euler, gammaRef, tol: 1e-15))
        #expect(approxEqual(Constants<Double>.catalan, catalanRef, tol: 1e-15))
        #expect(approxEqual(Constants<Double>.zetaThree, zeta3Ref, tol: 1e-15))
    }

    @Test
    func goldenRatio() {
        let phi = Constants<Double>.phi
        let constructed = (1.0 + (5.0 as Double).squareRoot()) / 2.0
        #expect(approxEqual(phi, constructed, tol: tol))
        #expect(phi > 1 && phi < 2)
    }
}

// MARK: - Float tests

@Suite("Constants<Float> tests")
struct FloatConstantsTests {
    let tol: Float = 1e-6

    @Test
    func piFamilyAndDerived() {
        let pi = Constants<Float>.pi
        #expect(pi.isFinite && pi > 3 && pi < 4)

        #expect(approxEqual(Constants<Float>.twoPi, 2 * pi, tol: tol))
        #expect(approxEqual(Constants<Float>.halfPi, pi / 2, tol: tol))
        #expect(approxEqual(Constants<Float>.quarterPi, pi / 4, tol: tol))
        #expect(approxEqual(Constants<Float>.thirdPi, pi / 3, tol: tol))
        #expect(approxEqual(Constants<Float>.twoThirdsPi, 2 * pi / 3, tol: tol))
        #expect(approxEqual(Constants<Float>.threeQuartersPi, 3 * pi / 4, tol: tol))
        #expect(approxEqual(Constants<Float>.sixthPi, pi / 6, tol: tol))
        #expect(approxEqual(Constants<Float>.piSqr, pi * pi, tol: tol))
    }

    @Test
    func rootsAndReciprocals() {
        let pi = Constants<Float>.pi
        let rootTwo = Constants<Float>.rootTwo
        let rootThree = Constants<Float>.rootThree
        let rootPi = Constants<Float>.rootPi

        #expect(approxEqual(rootTwo * rootTwo, 2.0, tol: tol))
        #expect(approxEqual(rootThree * rootThree, 3.0, tol: tol))
        #expect(approxEqual(rootPi * rootPi, pi, tol: tol))

        #expect(approxEqual(Constants<Float>.oneDivPi, 1.0 / pi, tol: tol))
        #expect(approxEqual(Constants<Float>.oneDivTwoPi, 1.0 / (2.0 * pi), tol: tol))
        #expect(approxEqual(Constants<Float>.oneDivRootPi, 1.0 / rootPi, tol: tol))
        #expect(approxEqual(Constants<Float>.twoDivPi, 2.0 / pi, tol: tol))
        #expect(approxEqual(Constants<Float>.twoDivRootPi, 2.0 / rootPi, tol: tol))
    }

    @Test
    func logsAndExponential() {
        let ln2 = Constants<Float>.lnTwo
        let ln10 = Constants<Float>.lnTen
        let lnln2 = Constants<Float>.lnLnTwo
        let e = Constants<Float>.e

        #expect(approxEqual(ln2, Darwin.log(2.0), tol: tol))
        #expect(approxEqual(ln10, Darwin.log(10.0), tol: tol))
        #expect(approxEqual(lnln2, Darwin.log(ln2), tol: tol))
        #expect(approxEqual(e, Darwin.exp(1.0), tol: tol))
    }

    @Test
    func specialConstants() {
        let gammaRef: Float = 0.5772156649015329
        let catalanRef: Float = 0.9159655941772190
        let zeta3Ref: Float = 1.2020569031595943

        #expect(approxEqual(Constants<Float>.euler, gammaRef, tol: tol))
        #expect(approxEqual(Constants<Float>.catalan, catalanRef, tol: tol))
        #expect(approxEqual(Constants<Float>.zetaThree, zeta3Ref, tol: tol))
    }

    @Test
    func goldenRatio() {
        let phi = Constants<Float>.phi
        let constructed = (1.0 as Float + (5.0 as Float).squareRoot()) / 2.0
        #expect(approxEqual(phi, constructed, tol: tol))
        #expect(phi > 1 && phi < 2)
    }
}

// MARK: - Float80 tests (x86_64 only)

#if arch(x86_64)
@Suite("Constants<Float80> tests")
struct Float80ConstantsTests {
    let tol: Float80 = 1e-18

    @Test
    func piFamilyAndDerived() {
        let pi = Constants<Float80>.pi
        #expect(pi.isFinite && pi > 3 && pi < 4)

        #expect(approxEqual(Constants<Float80>.twoPi, 2 * pi, tol: tol))
        #expect(approxEqual(Constants<Float80>.halfPi, pi / 2, tol: tol))
        #expect(approxEqual(Constants<Float80>.quarterPi, pi / 4, tol: tol))
        #expect(approxEqual(Constants<Float80>.thirdPi, pi / 3, tol: tol))
        #expect(approxEqual(Constants<Float80>.twoThirdsPi, 2 * pi / 3, tol: tol))
        #expect(approxEqual(Constants<Float80>.threeQuartersPi, 3 * pi / 4, tol: tol))
        #expect(approxEqual(Constants<Float80>.sixthPi, pi / 6, tol: tol))
        #expect(approxEqual(Constants<Float80>.piSqr, pi * pi, tol: tol))
    }

    @Test
    func rootsAndReciprocals() {
        let pi = Constants<Float80>.pi
        let rootTwo = Constants<Float80>.rootTwo
        let rootThree = Constants<Float80>.rootThree
        let rootPi = Constants<Float80>.rootPi

        #expect(approxEqual(rootTwo * rootTwo, 2.0, tol: tol))
        #expect(approxEqual(rootThree * rootThree, 3.0, tol: tol))
        #expect(approxEqual(rootPi * rootPi, pi, tol: tol))

        #expect(approxEqual(Constants<Float80>.oneDivPi, 1.0 / pi, tol: tol))
        #expect(approxEqual(Constants<Float80>.oneDivTwoPi, 1.0 / (2.0 * pi), tol: tol))
        #expect(approxEqual(Constants<Float80>.oneDivRootPi, 1.0 / rootPi, tol: tol))
        #expect(approxEqual(Constants<Float80>.twoDivPi, 2.0 / pi, tol: tol))
        #expect(approxEqual(Constants<Float80>.twoDivRootPi, 2.0 / rootPi, tol: tol))
    }

    @Test
    func logsAndExponential() {
        let ln2 = Constants<Float80>.lnTwo
        let ln10 = Constants<Float80>.lnTen
        let lnln2 = Constants<Float80>.lnLnTwo
        let e = Constants<Float80>.e

        #expect(approxEqual(ln2, .log(2.0), tol: tol))
        #expect(approxEqual(ln10, .log(10.0), tol: tol))
        #expect(approxEqual(lnln2, .log(ln2), tol: tol))
        #expect(approxEqual(e, .exp(1.0), tol: tol))
    }

    @Test
    func specialConstants() {
        let gammaRef: Float80 = 0.5772156649015328606
        let catalanRef: Float80 = 0.91596559417721901505
        let zeta3Ref: Float80 = 1.2020569031595942854

        #expect(approxEqual(Constants<Float80>.euler, gammaRef, tol: tol))
        #expect(approxEqual(Constants<Float80>.catalan, catalanRef, tol: tol))
        #expect(approxEqual(Constants<Float80>.zetaThree, zeta3Ref, tol: tol))
    }

    @Test
    func goldenRatio() {
        let phi = Constants<Float80>.phi
        let constructed = (1.0 as Float80 + (5.0 as Float80).squareRoot()) / 2.0
        #expect(approxEqual(phi, constructed, tol: tol))
        #expect(phi > 1 && phi < 2)
    }
}
#endif


