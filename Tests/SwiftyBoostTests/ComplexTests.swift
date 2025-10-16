//
//  Created by VT on 15.10.25.
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
import Foundation
@testable import SwiftyBoost

@Suite("Complex basic arithmetic and properties (Double)")
struct ComplexDoubleBasicTests {

    @Test("Construction and constants")
    func constructionAndConstants() {
        let z = Complex<Double>(3, -4)
        #expect(z.real == 3)
        #expect(z.imag == -4)
        #expect(Complex<Double>.zero == Complex(0, 0))
        #expect(Complex<Double>.one == Complex(1, 0))
        #expect(Complex<Double>.i == Complex(0, 1))
    }

    @Test("Addition, subtraction, negation")
    func addSubNeg() {
        let a = Complex(3.0, 2.0)
        let b = Complex(-1.0, 4.0)
        #expect(a + b == Complex(2.0, 6.0))
        #expect(a - b == Complex(4.0, -2.0))
        #expect(-a == Complex(-3.0, -2.0))
    }

    @Test("Multiplication")
    func multiply() {
        let a = Complex(3.0, 2.0)
        let b = Complex(-1.0, 4.0)
        // (3+2i)(-1+4i) = -3 + 12i - 2i + 8i^2 = (-11) + 10i
        #expect(a * b == Complex(-11.0, 10.0))
    }

    @Test("Division (Smith’s algorithm) and reciprocal")
    func divideAndReciprocal() {
        let a = Complex(3.0, 2.0)
        let b = Complex(-1.0, 4.0)

        let q = a / b
        // Validate by multiplying back
        let back = q * b
        #expect(abs(back.real - a.real) < 1e-12)
        #expect(abs(back.imag - a.imag) < 1e-12)

        // Reciprocal: z * (1/z) = 1
        let r = b.reciprocal()
        let one = b * r
        #expect(abs(one.real - 1.0) < 1e-12)
        #expect(abs(one.imag - 0.0) < 1e-12)
    }

    @Test("Mixed scalar arithmetic")
    func mixedScalar() {
        let z = Complex(3.0, -2.0)
        #expect(z + 2.0 == Complex(5.0, -2.0))
        #expect(2.0 + z == Complex(5.0, -2.0))
        #expect(z - 1.5 == Complex(1.5, -2.0))
        #expect(1.0 - z == Complex(-2.0, 2.0))
        #expect(z * 2.0 == Complex(6.0, -4.0))
        #expect(2.0 * z == Complex(6.0, -4.0))
        #expect(z / 2.0 == Complex(1.5, -1.0))
        // scalar / complex
        let q = 2.0 / z
        #expect(abs((q * z).real - 2.0) < 1e-12)
        #expect(abs((q * z).imag - 0.0) < 1e-12)
    }

    @Test("Magnitude, norm, conjugate, normalization")
    func magnitudeEtc() {
        let z = Complex(3.0, -4.0)
        #expect(z.norm == 25.0)
        #expect(abs(z.magnitude - 5.0) < 1e-12)
        #expect(z.conjugate == Complex(3.0, 4.0))
        let n = z.normalized()
        #expect(abs(n.magnitude - 1.0) < 1e-12)
        // Direction preserved (real/imag ratio)
        #expect(abs(n.real / n.imag - z.real / z.imag) < 1e-12)
    }

    @Test("Integer powers by squaring")
    func integerPowers() {
        let z = Complex(1.0, 1.0)
        #expect(z.pow(0) == Complex.one)
        #expect(z.pow(1) == z)
        // (1+i)^2 = 2i
        #expect(z.pow(2) == Complex(0.0, 2.0))
        // Negative power: z^(-1) * z = 1
        let inv = z.pow(-1)
        let back = inv * z
        #expect(abs(back.real - 1.0) < 1e-12)
        #expect(abs(back.imag - 0.0) < 1e-12)
    }

    @Test("Square root principal branch")
    func squareRootPrincipal() {
        // sqrt(4) = 2
        let a = Complex(4.0, 0.0).squareRoot
        #expect(abs(a.real - 2.0) < 1e-12)
        #expect(abs(a.imag - 0.0) < 1e-12)

        // sqrt(-4) = 2i
        let b = Complex(-4.0, 0.0).squareRoot
        #expect(abs(b.real - 0.0) < 1e-12)
        #expect(abs(b.imag - 2.0) < 1e-12)

        // General case: check squaring back
        let z = Complex(3.0, 4.0)
        let s = z.squareRoot
        let back = s * s
        #expect(abs(back.real - z.real) < 1e-10)
        #expect(abs(back.imag - z.imag) < 1e-10)
    }

    @Test("Polar construction and phase")
    func polarAndPhase() {
        let r: Double = 2.0
        let theta = Constants<Double>.pi / 4.0 // 45°
        let z = Complex<Double>.fromPolar(radius: r, phase: theta)
        // Expected: (sqrt(2), sqrt(2))
        let rt2 = r * sqrt(0.5)
        #expect(abs(z.real - rt2) < 1e-12)
        #expect(abs(z.imag - rt2) < 1e-12)
        #expect(abs(z.phase - theta) < 1e-12)
    }

    @Test("Exp and Log (Double specialization)")
    func expAndLog() {
        // exp(0) = 1
        let z0 = Complex(0.0, 0.0).exp
        #expect(abs(z0.real - 1.0) < 1e-12)
        #expect(abs(z0.imag - 0.0) < 1e-12)

        // exp(iπ) = -1
        let ziPi = Complex(0.0, .pi).exp
        #expect(abs(ziPi.real + 1.0) < 1e-12)
        #expect(abs(ziPi.imag - 0.0) < 1e-10) // allow small sin(pi) noise

        // log of positive real: log(4) = ln 4 + i 0
        let l = Complex(4.0, 0.0).log
        #expect(abs(l.real - log(4.0)) < 1e-12)
        #expect(abs(l.imag - 0.0) < 1e-12)
    }

    @Test("Codable round-trip JSON")
    func codableRoundTrip() throws {
        let z = Complex(3.25, -4.75)
        let data = try JSONEncoder().encode(z)
        let decoded = try JSONDecoder().decode(Complex<Double>.self, from: data)
        #expect(decoded == z)
    }

    @Test("Hashable semantics")
    func hashable() {
        let a = Complex(1.0, 2.0)
        let b = Complex(1.0, 2.0)
        let c = Complex(2.0, 1.0)
        var set = Set<Complex<Double>>()
        set.insert(a)
        set.insert(b)
        set.insert(c)
        #expect(set.count == 2)
        #expect(set.contains(a))
        #expect(set.contains(c))
    }

    @Test("Description formatting")
    func description() {
        let z1 = Complex(3.0, 2.0)
        #expect(z1.description == "3.0 + 2.0i" || z1.description == "3 + 2i")
        let z2 = Complex(3.0, -2.0)
        #expect(z2.description == "3.0 - 2.0i" || z2.description == "3 - 2i")
    }
}

@Suite("Complex Float specializations")
struct ComplexFloatTests {

    @Test("Phase and fromPolar (Float)")
    func phaseAndPolarFloat() {
        let r: Float = 2
        let theta: Float = .pi / 6 // 30°
        let z = Complex<Float>.fromPolar(radius: r, phase: theta)
        #expect(abs(z.phase - theta) < 1e-5)

        // Check back to radius via magnitude
        #expect(abs(z.magnitude - r) < 1e-5)
    }

    @Test("Exp and Log (Float specialization)")
    func expLogFloat() {
        let z0 = Complex<Float>(0, 0).exp
        #expect(abs(z0.real - 1) < 1e-6)
        #expect(abs(z0.imag - 0) < 1e-6)

        let ziPi = Complex<Float>(0, .pi).exp
        #expect(abs(ziPi.real + 1) < 1e-5)
        #expect(abs(ziPi.imag - 0) < 1e-4)

        let l = Complex<Float>(4, 0).log
        #expect(abs(l.real - logf(4)) < 1e-6)
        #expect(abs(l.imag - 0) < 1e-6)
    }
}

#if arch(x86_64)
@Suite("Complex Float80 specializations (x86_64)")
struct ComplexFloat80Tests {

    @Test("Phase and fromPolar (Float80)")
    func phaseAndPolarFloat80() {
        let r: Float80 = 3
        let theta: Float80 = .pi / 3
        let z = Complex<Float80>.fromPolar(radius: r, phase: theta)
        #expect(abs(z.phase - theta) < 1e-12)
        #expect(abs(z.magnitude - r) < 1e-12)
    }

    @Test("Exp and Log (Float80 specialization)")
    func expLogFloat80() {
        let z0 = Complex<Float80>(0, 0).exp
        #expect(abs(z0.real - 1) < 1e-12)
        #expect(abs(z0.imag - 0) < 1e-12)

        let ziPi = Complex<Float80>(0, .pi).exp
        #expect(abs(ziPi.real + 1) < 1e-12)
        #expect(abs(ziPi.imag - 0) < 1e-10)

        let l = Complex<Float80>(4, 0).log
        #expect(abs(l.real - logl(4)) < 1e-12)
        #expect(abs(l.imag - 0) < 1e-12)
    }
}
#endif
