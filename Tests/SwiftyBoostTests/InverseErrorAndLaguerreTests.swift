// Tests/SwiftyBoostTests/InverseErrorAndLaguerreTests.swift
import Testing
import SwiftyBoost

@Suite("Inverse error function tests (Double/Float/Float80 where available)")
struct InverseErrorFunctionTests {
    // Representative probabilities in (-1, 1)
    let psD: [Double] = [-0.999999, -0.99, -0.9, -0.5, -0.1, 0.0, 0.1, 0.5, 0.9, 0.99, 0.999999]
    let psF: [Float]  = [-0.99, -0.9, -0.5, -0.1, 0.0, 0.1, 0.5, 0.9, 0.99]

    @Test("Round-trip erf(erfInv(p)) ≈ p (Double)")
    func roundTripDouble() async throws {
        for p in psD {
            let x: Double = try SpecialFunctions.inverseErrorFunction(p)
            let r: Double = try SpecialFunctions.errorFunction(x)
            #expect(abs(r - p) <= 2e-12, "Round-trip failed for p=\(p), got \(r)")
        }
    }

    @Test("Round-trip erf(erfInv(p)) ≈ p (Float)")
    func roundTripFloat() async throws {
        for p in psF {
            let x: Float = try SpecialFunctions.inverseErrorFunction(p)
            let r: Float = try SpecialFunctions.errorFunction(x)
            #expect(abs(r - p) <= 1e-6, "Round-trip failed for p=\(p), got \(r)")
        }
    }

    #if arch(x86_64)
    @Test("Round-trip erf(erfInv(p)) ≈ p (Float80)")
    func roundTripFloat80() async throws {
        let psL: [Float80] = [-0.999999, -0.99, -0.9, -0.5, -0.1, 0.0, 0.1, 0.5, 0.9, 0.99, 0.999999, Float80(1.5).nextUp]
        for p in psL {
            let x: Float80 = try SpecialFunctions.inverseErrorFunction(p)
            let r: Float80 = try SpecialFunctions.errorFunction(x)
            #expect(abs(r - p) <= 5e-18, "Round-trip failed for p=\(p), got \(r)")
        }
    }
    #endif

    @Test("erfc(erfInv(p)) ≈ 1 - p (Double)")
    func erfcRelationDouble() async throws {
        for p in psD {
            let x: Double = try SpecialFunctions.inverseErrorFunction(p)
            let r: Double = try SpecialFunctions.complementaryErrorFunction(x)
            #expect(abs(r - (1 - p)) <= 2e-12, "erfc relation failed for p=\(p), got \(r)")
        }
    }

    @Test("Domain errors for inverseErrorFunction (Double)")
    func domainErrorsDouble() async throws {
        // Non-finite
        #expect(throws: SpecialFunctionError.parameterNotFinite(name: "z")) {
            _ = try SpecialFunctions.inverseErrorFunction(Double.nan)
        }
        #expect(throws: SpecialFunctionError.parameterNotFinite(name: "z")) {
            _ = try SpecialFunctions.inverseErrorFunction(Double.infinity)
        }
        #expect(throws: SpecialFunctionError.parameterNotFinite(name: "z")) {
            _ = try SpecialFunctions.inverseErrorFunction(-Double.infinity)
        }
        // |z| >= 1
        #expect(throws: SpecialFunctionError.parameterOutOfRange(name: "z", min: (-1.0).nextUp, max: (1.0).nextDown)) {
            _ = try SpecialFunctions.inverseErrorFunction(1.0)
        }
        #expect(throws: SpecialFunctionError.parameterOutOfRange(name: "z", min: (-1.0).nextUp, max: (1.0).nextDown)) {
            _ = try SpecialFunctions.inverseErrorFunction(-1.0)
        }
        #expect(throws: SpecialFunctionError.parameterOutOfRange(name: "z", min: (-1.0).nextUp, max: (1.0).nextDown)) {
            _ = try SpecialFunctions.inverseErrorFunction(1.0000001)
        }
        #expect(throws: SpecialFunctionError.parameterOutOfRange(name: "z", min: (-1.0).nextUp, max: (1.0).nextDown)) {
            _ = try SpecialFunctions.inverseErrorFunction(-1.0000001)
        }
    }

    @Test("Domain errors for inverseErrorFunction (Float)")
    func domainErrorsFloat() async throws {
        // Non-finite
        #expect(throws: SpecialFunctionError.parameterNotFinite(name: "z")) {
            _ = try SpecialFunctions.inverseErrorFunction(Float.nan)
        }
        #expect(throws: SpecialFunctionError.parameterNotFinite(name: "z")) {
            _ = try SpecialFunctions.inverseErrorFunction(Float.infinity)
        }
        #expect(throws: SpecialFunctionError.parameterNotFinite(name: "z")) {
            _ = try SpecialFunctions.inverseErrorFunction(-Float.infinity)
        }
        // |z| >= 1
        #expect(throws: SpecialFunctionError.parameterOutOfRange(name: "z", min: Double(Float(-1).nextUp), max: Double(Float(1).nextDown))) {
            _ = try SpecialFunctions.inverseErrorFunction(Float(1))
        }
        #expect(throws: SpecialFunctionError.parameterOutOfRange(name: "z", min: Double(Float(-1).nextUp), max: Double(Float(1).nextDown))) {
            _ = try SpecialFunctions.inverseErrorFunction(Float(-1))
        }
    }
}

@Suite("Laguerre polynomial tests (Double/Float/Float80 where available)")
struct LaguerreTests {

    // Closed-form small-n Laguerre polynomials:
    // L0(x) = 1
    // L1(x) = 1 - x
    // L2(x) = 1 - 2x + x^2/2
    // L3(x) = 1 - 3x + 3x^2/2 - x^3/6

    @Test("L_n(x) closed forms for n = 0...3 (Double)")
    func basicLaguerreDouble() async throws {
        let xs: [Double] = [0, 0.5, 1, 2, -1]
        for x in xs {
            #expect(try SpecialFunctions.laguerre(0, x) == 1)
            #expect(abs(try SpecialFunctions.laguerre(1, x) - (1 - x)) <= 1e-15)
            #expect(abs(try SpecialFunctions.laguerre(2, x) - (1 - 2*x + 0.5*x*x)) <= 2e-15)
            let l3exp = 1 - 3*x + 1.5*x*x - (1.0/6.0)*x*x*x
            #expect(abs(try SpecialFunctions.laguerre(3, x) - l3exp) <= 8e-15)
        }
    }

    @Test("L_n(x) closed forms for n = 0...3 (Float)")
    func basicLaguerreFloat() async throws {
        let xs: [Float] = [0, 0.5, 1, 2, -1]
        for x in xs {
            #expect(try SpecialFunctions.laguerre(0, x) == 1)
            #expect(abs(try SpecialFunctions.laguerre(1, x) - (1 - x)) <= 2e-6)
            #expect(abs(try SpecialFunctions.laguerre(2, x) - (1 - 2*x + 0.5*x*x)) <= 3e-6)
            let l3exp = 1 - 3*x + 1.5*x*x - (1.0/6.0)*x*x*x
            #expect(abs(try SpecialFunctions.laguerre(3, x) - l3exp) <= 6e-6)
        }
    }

    #if arch(x86_64)
    @Test("L_n(x) closed forms for n = 0...3 (Float80)")
    func basicLaguerreFloat80() async throws {
        let xs: [Float80] = [0, 0.5, 1, 2, -1]
        for x in xs {
            #expect(try SpecialFunctions.laguerre(0, x) == 1)
            #expect(abs(try SpecialFunctions.laguerre(1, x) - (1 - x)) <= 1e-18)
            #expect(abs(try SpecialFunctions.laguerre(2, x) - (1 - 2*x + 0.5*x*x)) <= 2e-18)
            let l3exp: Float80 = 1 - 3*x + 1.5*x*x - (1.0/6.0)*x*x*x
            #expect(abs(try SpecialFunctions.laguerre(3, x) - l3exp) <= 6e-18)
        }
    }
    #endif

    @Test("Associated Laguerre: L_n^0(x) = L_n(x) and L_n^m(0) = C(n+m, n) (Double)")
    func assocLaguerreIdentitiesDouble() async throws {
        // Check L_n^0(x) = L_n(x)
        let xs: [Double] = [-1, 0, 0.25, 1.25, 3]
        for n in 0...6 {
            for x in xs {
                let a = try SpecialFunctions.assocLaguerre(n, 0, x)
                let b = try SpecialFunctions.laguerre(n, x)
                #expect(abs(a - b) <= 1e-12, "Mismatch at n=\(n), x=\(x)")
            }
        }
        // Check L_n^m(0) = binomial(n+m, n)
        for n in 0...8 {
            for m in 0...8 {
                let lhs = try SpecialFunctions.assocLaguerre(n, m, 0.0 as Double)
                let rhs: Double = try SpecialFunctions.binomial_coeff(UInt32(n + m), UInt32(n))
                #expect(abs(lhs - rhs) <= 1e-12, "Mismatch at n=\(n), m=\(m)")
            }
        }
    }

    @Test("Associated Laguerre: L_n^0(x) = L_n(x) and L_n^m(0) binomial (Float)")
    func assocLaguerreIdentitiesFloat() async throws {
        let xs: [Float] = [-1, 0, 0.25, 1.25, 3]
        for n in 0...6 {
            for x in xs {
                let a = try SpecialFunctions.assocLaguerre(n, 0, x)
                let b = try SpecialFunctions.laguerre(n, x)
                #expect(abs(a - b) <= 2e-5, "Mismatch at n=\(n), x=\(x)")
            }
        }
        for n in 0...8 {
            for m in 0...8 {
                let lhs = try SpecialFunctions.assocLaguerre(n, m, 0.0 as Float)
                let rhs: Float = try SpecialFunctions.binomial_coeff_f(UInt32(n + m), UInt32(n))
                #expect(abs(lhs - rhs) <= 2e-5, "Mismatch at n=\(n), m=\(m)")
            }
        }
    }

    #if arch(x86_64)
    @Test("Associated Laguerre: identities (Float80)")
    func assocLaguerreIdentitiesFloat80() async throws {
        let xs: [Float80] = [-1, 0, 0.25, 1.25, 3]
        for n in 0...6 {
            for x in xs {
                let a = try SpecialFunctions.assocLaguerre(n, 0, x)
                let b = try SpecialFunctions.laguerre(n, x)
                #expect(abs(a - b) <= 1e-18, "Mismatch at n=\(n), x=\(x)")
            }
        }
        for n in 0...8 {
            for m in 0...8 {
                let lhs = try SpecialFunctions.assocLaguerre(n, m, 0.0 as Float80)
                let rhs: Float80 = try SpecialFunctions.binomial_coeff_l(UInt32(n + m), UInt32(n))
                #expect(abs(lhs - rhs) <= 1e-18, "Mismatch at n=\(n), m=\(m)")
            }
        }
    }
    #endif

    @Test("Laguerre error paths (Double)")
    func laguerreErrorsDouble() async throws {
        // n < 0
        #expect(throws: SpecialFunctionError.parameterNotPositive(name: "n")) {
            _ = try SpecialFunctions.laguerre(-1, 0.5 as Double)
        }
        // m < 0
        #expect(throws: SpecialFunctionError.parameterNotPositive(name: "m")) {
            _ = try SpecialFunctions.assocLaguerre(2, -3, 1.0 as Double)
        }
        // non-finite x
        #expect(throws: SpecialFunctionError.parameterNotFinite(name: "x")) {
            _ = try SpecialFunctions.laguerre(1, Double.nan)
        }
        #expect(throws: SpecialFunctionError.parameterNotFinite(name: "x")) {
            _ = try SpecialFunctions.assocLaguerre(1, 1, Double.infinity)
        }
    }

    @Test("Laguerre error paths (Float)")
    func laguerreErrorsFloat() async throws {
        #expect(throws: SpecialFunctionError.parameterNotPositive(name: "n")) {
            _ = try SpecialFunctions.laguerre(-1, 0.5 as Float)
        }
        #expect(throws: SpecialFunctionError.parameterNotPositive(name: "m")) {
            _ = try SpecialFunctions.assocLaguerre(2, -3, 1.0 as Float)
        }
        #expect(throws: SpecialFunctionError.parameterNotFinite(name: "x")) {
            _ = try SpecialFunctions.laguerre(1, Float.nan)
        }
        #expect(throws: SpecialFunctionError.parameterNotFinite(name: "x")) {
            _ = try SpecialFunctions.assocLaguerre(1, 1, Float.infinity)
        }
    }

    #if arch(x86_64)
    @Test("Laguerre error paths (Float80)")
    func laguerreErrorsFloat80() async throws {
        #expect(throws: SpecialFunctionError.parameterNotPositive(name: "n")) {
            _ = try SpecialFunctions.laguerre(-1, 0.5 as Float80)
        }
        #expect(throws: SpecialFunctionError.parameterNotPositive(name: "m")) {
            _ = try SpecialFunctions.assocLaguerre(2, -3, 1.0 as Float80)
        }
        #expect(throws: SpecialFunctionError.parameterNotFinite(name: "x")) {
            _ = try SpecialFunctions.laguerre(1, Float80.nan)
        }
        #expect(throws: SpecialFunctionError.parameterNotFinite(name: "x")) {
            _ = try SpecialFunctions.assocLaguerre(1, 1, Float80.infinity)
        }
    }
    #endif
}
