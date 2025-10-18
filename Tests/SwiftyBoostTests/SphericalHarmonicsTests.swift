//
//  Created by VT on 16.10.25.
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
//
//  SphericalHarmonicsTests.swift
//  SwiftyBoostTests
//
//  Created by VT on 16.10.25.
//

import Testing
@testable import SwiftyBoost

@Suite("Spherical Harmonics (Double)")
struct SphericalHarmonicsDoubleTests {

    // Helper: relative tolerance check
    @inlinable
    func close(
        _ a: Double,
        _ b: Double,
        rel tol: Double = 1e-12,
        fileID: String = #fileID,
        filePath: String = #filePath,
        line: Int = #line,
        column: Int = #column
    ) {
        let da = abs(a - b)
        let scale = max(1.0, abs(a), abs(b))
        #expect(
            da <= tol * scale,
            "Expected \(a) ≈ \(b) within rel tol \(tol), got |Δ|=\(da)",
            sourceLocation: .init(fileID: fileID, filePath: filePath, line: line, column: column)
        )
    }

    @Test("Spot checks against Boost reference values (real part)")
    func spotReal() throws {
        // From Boost's test_spots:
        // spherical_harmonic_r(3, 2, 0.5, 0) ≈ 0.20614605996878713306922867918
        close(
            SpecialFunctions.Y_real(n: 3, m: 2, theta: 0.5, phi: 0.0),
            0.2061460599687871330692286791802688341213
        )

        // spherical_harmonic_r(20, 10, 0.75, -0.25) ≈ 0.06197787102219208244
        close(
            SpecialFunctions.Y_real(n: 20, m: 10, theta: 0.75, phi: -0.25),
            0.06197787102219208244041677775577045124092
        )

        // spherical_harmonic_r(40, 15, -0.75, 2.25) ≈ 0.28069048250457456873
        close(
            SpecialFunctions.Y_real(n: 40, m: 15, theta: -0.75, phi: 2.25),
            0.2806904825045745687343492963236868973484
        )

        // spherical_harmonic_r(20, 14, -0.75, 2.25) ≈ 0.34792181861334354667
        close(
            SpecialFunctions.Y_real(n: 20, m: 14, theta: -0.75, phi: 2.25),
            0.3479218186133435466692822481919867452442
        )
    }

    @Test("Spot checks against Boost reference values (imag part)")
    func spotImag() throws {
        // spherical_harmonic_i(20, 10, 0.75, -0.25) ≈ 0.04629885158895932341
        close(
            SpecialFunctions.Y_imag(n: 20, m: 10, theta: 0.75, phi: -0.25),
            0.04629885158895932341185988759669916977920
        )

        // spherical_harmonic_i(40, 15, -0.75, 2.25) ≈ -0.29339184446566035823
        close(
            SpecialFunctions.Y_imag(n: 40, m: 15, theta: -0.75, phi: 2.25),
            -0.2933918444656603582282372590387544902135
        )

        // spherical_harmonic_i(20, 14, -0.75, 2.25) ≈ 0.02932010666852638796
        close(
            SpecialFunctions.Y_imag(n: 20, m: 14, theta: -0.75, phi: 2.25),
            0.0293201066685263879566422194539567289974
        )
    }

    @Test("Negative m identity: Y_l^{-m} = (-1)^m * conj(Y_l^{m})")
    func negativeMIdentity() throws {
        let l: UInt = 40
        let m = 15
        let theta = -0.75
        let phi = 2.25

        let yp = SpecialFunctions.sphericalHarmonic(n: l, m: m, theta: theta, phi: phi)
        let yn = SpecialFunctions.sphericalHarmonic(n: l, m: -m, theta: theta, phi: phi)

        // Branchless computation of (-1)^m to avoid "will never be executed" warning.
        let sign = 1.0 - 2.0 * Double(m & 1)
        let conjYp = ComplexD(yp.real, -yp.imag)
        let rhs = ComplexD(sign * conjYp.real, sign * conjYp.imag)

        close(yn.real, rhs.real)
        close(yn.imag, rhs.imag)
    }
}

@Suite("Spherical Harmonics (Float)")
struct SphericalHarmonicsFloatTests {

    @inlinable
    func close(
        _ a: Float,
        _ b: Float,
        rel tol: Float = 1e-5,
        fileID: String = #fileID,
        filePath: String = #filePath,
        line: Int = #line,
        column: Int = #column
    ) {
        let da = abs(a - b)
        let scale = max(Float(1), abs(a), abs(b))
        #expect(
            da <= tol * scale,
            "Expected \(a) ≈ \(b) within rel tol \(tol), got |Δ|=\(da)",
            sourceLocation: .init(fileID: fileID, filePath: filePath, line: line, column: column)
        )
    }

    @Test("Spot checks (Float, via Double path)")
    func spotFloat() throws {
        // Reuse a couple of double cases, cast to Float, allow looser tolerance.
        let r1 = SpecialFunctions.Y(n: 3, m: 2, theta: Float(0.5), phi: Float(0.0))
        close(r1.real, Float(0.20614605996878713))
        close(r1.imag, 0.0)

        let r2 = SpecialFunctions.Y(n: 20, m: 10, theta: Float(0.75), phi: Float(-0.25))
        close(r2.real, Float(0.06197787102219208))
        close(r2.imag, Float(0.04629885158895932))
    }
}

#if arch(x86_64)
@Suite("Spherical Harmonics (Float80)")
struct SphericalHarmonicsFloat80Tests {

    @inlinable
    func close(
        _ a: Float80,
        _ b: Float80,
        rel tol: Float80 = 1e-12,
        fileID: String = #fileID,
        filePath: String = #filePath,
        line: Int = #line,
        column: Int = #column
    ) {
        let diff = a - b
        let da: Float80 = diff < 0 ? -diff : diff
        let aa: Float80 = a < 0 ? -a : a
        let bb: Float80 = b < 0 ? -b : b
        let scale = max(Float80(1), aa, bb)
        #expect(
            da <= tol * scale,
            "Expected \(a) ≈ \(b) within rel tol \(tol), got |Δ|=\(da)",
            sourceLocation: .init(fileID: fileID, filePath: filePath, line: line, column: column)
        )
    }

    @Test("Spot checks (Float80, via Double path)")
    func spotFloat80() throws {
        let r = SpecialFunctions.Y(n: 40, m: 15, theta: Float80(-0.75), phi: Float80(2.25))
        close(r.real, 0.2806904825045745687343)
        close(r.imag, -0.2933918444656603582282)
    }
}
#endif
