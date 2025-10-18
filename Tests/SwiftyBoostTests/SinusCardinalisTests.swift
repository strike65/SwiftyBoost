import Testing
@testable import SwiftyBoost
import Foundation

@Suite("Sinus cardinalis (sinc/sinhc) tests")
struct SinusCardinalisTests {

    // MARK: - sinc_pi (Double)
    @Test("sinc_pi basic identities (Double)")
    func sincPiDouble() throws {
        // x = 0 → 1
        do {
            let v = SpecialFunctions.sinc_pi(0.0)
            #expect(abs(v - 1.0) <= 1e-15)
        }
        // x = 1 → sin(π)/π = 0
        do {
            let v = SpecialFunctions.sinc_pi(1.0)
            #expect(abs(v) <= 1e-15)
        }
        // x = 1/2 → 2/π
        let got = SpecialFunctions.sinc_pi(0.5)
        let exp = 2.0 / Double.pi
        #expect(abs(got - exp) <= 1e-15)
        // Even symmetry: sinc(−x) = sinc(x)
        let neg = SpecialFunctions.sinc_pi(-0.5)
        #expect(abs(neg - exp) <= 1e-15)
    }

    // MARK: - sinc_pi (Float)
    @Test("sinc_pi basic identities (Float)")
    func sincPiFloat() throws {
        do {
            let v0 = SpecialFunctions.sinc_pi(0 as Float)
            #expect(abs(v0 - 1) <= 1e-6)
        }
        do {
            let v1 = SpecialFunctions.sinc_pi(1 as Float)
            #expect(abs(v1) <= 1e-6)
        }
        let got = SpecialFunctions.sinc_pi(0.5 as Float)
        let exp: Float = 2.0 / Float.pi
        #expect(abs(got - exp) <= 1e-6)
        let neg = SpecialFunctions.sinc_pi(-0.5 as Float)
        #expect(abs(neg - exp) <= 1e-6)
    }

    // MARK: - sinc_pi (Float80)
    #if arch(x86_64)
    @Test("sinc_pi basic identities (Float80)")
    func sincPiFloat80() throws {
        do {
            let v0 = SpecialFunctions.sinc_pi(0 as Float80)
            #expect(abs(v0 - 1) <= 1e-16)
        }
        do {
            let v1 = SpecialFunctions.sinc_pi(1 as Float80)
            #expect(abs(v1) <= 1e-16)
        }
        let got = SpecialFunctions.sinc_pi(0.5 as Float80)
        let exp: Float80 = 2.0 / Float80.pi
        #expect(abs(got - exp) <= 1e-16)
        let neg = SpecialFunctions.sinc_pi(-0.5 as Float80)
        #expect(abs(neg - exp) <= 1e-16)
    }
    #endif

    // MARK: - sinhc_pi (Double)
    @Test("sinhc_pi basic identities (Double)")
    func sinhcPiDouble() throws {
        // x = 0 → 1 (limit)
        #expect(SpecialFunctions.sinhc_pi(0.0) == 1.0)
        // Compare to definition for representative points
        for x in [0.25, 0.5, 1.0] as [Double] {
            let got = SpecialFunctions.sinhc_pi(x)
            let exp = sinh(Double.pi * x) / (Double.pi * x)
            #expect(abs(got - exp) <= 1e-15, "Mismatch at x=\(x)")
        }
    }

    // MARK: - sinhc_pi (Float)
    @Test("sinhc_pi basic identities (Float)")
    func sinhcPiFloat() throws {
        #expect(SpecialFunctions.sinhc_pi(0 as Float) == 1)
        for x in [Float(0.25), 0.5, 1.0] {
            let got = SpecialFunctions.sinhc_pi(x)
            let exp: Float = sinh(Float.pi * x) / (Float.pi * x)
            #expect(abs(got - exp) <= 1e-6, "Mismatch at x=\(x)")
        }
    }

    // MARK: - sinhc_pi (Float80)
    #if arch(x86_64)
    @Test("sinhc_pi basic identities (Float80)")
    func sinhcPiFloat80() throws {
        #expect(SpecialFunctions.sinhc_pi(0 as Float80) == 1)
        for x in [Float80(0.25), 0.5, 1.0] {
            let got = SpecialFunctions.sinhc_pi(x)
            let exp: Float80 = sinh(Float80.pi * x) / (Float80.pi * x)
            #expect(abs(got - exp) <= 1e-16, "Mismatch at x=\(x)")
        }
    }
    #endif

    // MARK: - Complex relation: sinc(π i a) = sinhc(π a)
    @Test("sincc_pi(i a) == sinhc_pi(a) (Double)")
    func complexSincRelation() throws {
        let a: Double = 0.7
        let z = ComplexD(0.0, a)
        let s = SpecialFunctions.sincc_pi(z)
        let exp = SpecialFunctions.sinhc_pi(a)
        #expect(abs(s.imag) <= 1e-15, "Expected purely real value")
        #expect(abs(s.real - exp) <= 1e-14, "Mismatch for a=\(a)")
    }

    // MARK: - Complex: sincc_pi(z) = sin(π z)/(π z)
    @Test("sincc_pi(z) matches sin(pi z)/(pi z) (Double)")
    func complexSincMatchesClosedFormDouble() {
        let zs: [ComplexD] = [
            ComplexD(0.25, 0.0),
            ComplexD(0.3, 0.4),
            ComplexD(-0.6, 0.2)
        ]
        for z in zs {
            let pz = z * Double.pi
            let expected = pz.sin / pz
            let got = SpecialFunctions.sincc_pi(z)
            #expect(abs(got.real - expected.real) <= 1e-13, "Real mismatch for z=\(z)")
            #expect(abs(got.imag - expected.imag) <= 1e-13, "Imag mismatch for z=\(z)")
        }
    }

    @Test("sincc_pi(z) matches sin(pi z)/(pi z) (Float)")
    func complexSincMatchesClosedFormFloat() {
        let zs: [ComplexF] = [
            ComplexF(0.25, 0.0),
            ComplexF(0.3, 0.2),
            ComplexF(-0.5, 0.4)
        ]
        for z in zs {
            let pz = z * Float.pi
            let expected = pz.sin / pz
            let got = SpecialFunctions.sincc_pi(z)
            #expect(abs(got.real - expected.real) <= 1e-5, "Real mismatch for z=\(z)")
            #expect(abs(got.imag - expected.imag) <= 1e-5, "Imag mismatch for z=\(z)")
        }
    }

    #if arch(x86_64)
    @Test("sincc_pi(z) matches sin(pi z)/(pi z) (Float80)")
    func complexSincMatchesClosedFormFloat80() {
        let zs: [ComplexX] = [
            ComplexX(0.25, 0.0),
            ComplexX(0.3, 0.4),
            ComplexX(-0.6, 0.2)
        ]
        for z in zs {
            let pz = z * Float80.pi
            let expected = pz.sin / pz
            let got = SpecialFunctions.sincc_pi(z)
            #expect(abs(got.real - expected.real) <= 1e-16, "Real mismatch for z=\(z)")
            #expect(abs(got.imag - expected.imag) <= 1e-16, "Imag mismatch for z=\(z)")
        }
    }
    #endif
}
