// Tests/SwiftyBoostTests/ChebyshevTests.swift
import Testing
import SwiftyBoost
import Foundation

@Suite("Chebyshev polynomials (Double/Float/Float80 where available)")
struct ChebyshevTests {

    // MARK: - Closed forms and identities

    @Test("T_n(x) closed forms for n = 0...3 (Double)")
    func chebyshevTClosedFormsDouble() async throws {
        let xs: [Double] = [-1.5, -1, -0.5, 0, 0.3, 0.5, 1, 2]
        for x in xs {
            #expect(try SpecialFunctions.chebyshevT(0, x) == 1)
            #expect(abs(try SpecialFunctions.chebyshevT(1, x) - x) <= 1e-15)
            #expect(abs(try SpecialFunctions.chebyshevT(2, x) - (2*x*x - 1)) <= 3e-15)
            #expect(abs(try SpecialFunctions.chebyshevT(3, x) - (4*x*x*x - 3*x)) <= 2e-12)
        }
    }

    @Test("U_n(x) closed forms for n = 0...3 (Double)")
    func chebyshevUClosedFormsDouble() async throws {
        let xs: [Double] = [-1.5, -1, -0.5, 0, 0.3, 0.5, 1, 2]
        for x in xs {
            #expect(try SpecialFunctions.chebyshevU(0, x) == 1)
            #expect(abs(try SpecialFunctions.chebyshevU(1, x) - (2*x)) <= 2e-15)
            #expect(abs(try SpecialFunctions.chebyshevU(2, x) - (4*x*x - 1)) <= 4e-15)
            #expect(abs(try SpecialFunctions.chebyshevU(3, x) - (8*x*x*x - 4*x)) <= 1e-14)
        }
    }

    @Test("Trig identities on |x| ≤ 1 (Double)")
    func trigIdentitiesDouble() async throws {
        // Avoid θ near 0 or π to prevent sin θ ~ 0 for U_n identity.
        let thetas: [Double] = stride(from: 0.2, through: 2.8, by: 0.3).map { $0 } // radians
        let ns = [0, 1, 2, 3, 5, 7]
        for θ in thetas {
            let x = cos(θ)
            for n in ns {
                let Tn = try SpecialFunctions.chebyshevT(n, x)
                let Un = try SpecialFunctions.chebyshevU(n, x)
                let Tref = cos(Double(n) * θ)
                let Uref = sin(Double(n + 1) * θ) / sin(θ)
                #expect(abs(Tn - Tref) <= 2e-14, "T\(n)(x) trig identity failed at θ=\(θ)")
                #expect(abs(Un - Uref) <= 5e-14, "U\(n)(x) trig identity failed at θ=\(θ)")
            }
        }
    }

    @Test("Endpoint identities (Double)")
    func endpointsDouble() async throws {
        let ns = 0...12
        for n in ns {
            let T1 = try SpecialFunctions.chebyshevT(n, 1.0 as Double)
            let Tm1 = try SpecialFunctions.chebyshevT(n, -1.0 as Double)
            let U1 = try SpecialFunctions.chebyshevU(n, 1.0 as Double)
            let Um1 = try SpecialFunctions.chebyshevU(n, -1.0 as Double)
            let s = (n % 2 == 0) ? 1.0 : -1.0
            #expect(abs(T1 - 1.0) <= 1e-15)
            #expect(abs(Tm1 - s) <= 1e-15)
            #expect(abs(U1 - Double(n + 1)) <= 1e-14)
            #expect(abs(Um1 - s * Double(n + 1)) <= 1e-14)
        }
    }

    @Test("T_n(x) and U_n(x) closed forms (Float)")
    func chebyshevClosedFormsFloat() async throws {
        let xs: [Float] = [-1.5, -1, -0.5, 0, 0.3, 0.5, 1, 2]
        for x in xs {
            #expect(try SpecialFunctions.chebyshevT(0, x) == 1)
            #expect(abs(try SpecialFunctions.chebyshevT(1, x) - x) <= 2e-6)
            #expect(abs(try SpecialFunctions.chebyshevT(2, x) - (2*x*x - 1)) <= 4e-6)
            #expect(abs(try SpecialFunctions.chebyshevT(3, x) - (4*x*x*x - 3*x)) <= 8e-6)

            #expect(try SpecialFunctions.chebyshevU(0, x) == 1)
            #expect(abs(try SpecialFunctions.chebyshevU(1, x) - (2*x)) <= 3e-6)
            #expect(abs(try SpecialFunctions.chebyshevU(2, x) - (4*x*x - 1)) <= 6e-6)
            #expect(abs(try SpecialFunctions.chebyshevU(3, x) - (8*x*x*x - 4*x)) <= 1.2e-5)
        }
    }

    @Test("Trig identities on |x| ≤ 1 (Float)")
    func trigIdentitiesFloat() async throws {
        let thetas: [Float] = stride(from: Float(0.2), through: Float(2.8), by: Float(0.3)).map { $0 }
        let ns = [0, 1, 2, 3, 5, 7]
        for θ in thetas {
            let x = cos(θ)
            for n in ns {
                let Tn = try SpecialFunctions.chebyshevT(n, x)
                let Un = try SpecialFunctions.chebyshevU(n, x)
                let Tref = cos(Float(n) * θ)
                let Uref = sin(Float(n + 1) * θ) / sin(θ)
                #expect(abs(Tn - Tref) <= 1.5e-5)
                #expect(abs(Un - Uref) <= 2.5e-5)
            }
        }
    }

    @Test("Endpoint identities (Float)")
    func endpointsFloat() async throws {
        for n in 0...12 {
            let T1 = try SpecialFunctions.chebyshevT(n, 1.0 as Float)
            let Tm1 = try SpecialFunctions.chebyshevT(n, -1.0 as Float)
            let U1 = try SpecialFunctions.chebyshevU(n, 1.0 as Float)
            let Um1 = try SpecialFunctions.chebyshevU(n, -1.0 as Float)
            let s: Float = (n % 2 == 0) ? 1 : -1
            #expect(abs(T1 - 1) <= 2e-6)
            #expect(abs(Tm1 - s) <= 2e-6)
            #expect(abs(U1 - Float(n + 1)) <= 2e-5)
            #expect(abs(Um1 - s * Float(n + 1)) <= 2.5e-5)
        }
    }

    #if arch(x86_64)
    @Test("T_n(x) and U_n(x) closed forms (Float80)")
    func chebyshevClosedFormsFloat80() async throws {
        let xs: [Float80] = [-1.5, -1, -0.5, 0, 0.3, 0.5, 1, 2]
        for x in xs {
            #expect(try SpecialFunctions.chebyshevT(0, x) == 1)
            #expect(abs(try SpecialFunctions.chebyshevT(1, x) - x) <= 1e-18)
            #expect(abs(try SpecialFunctions.chebyshevT(2, x) - (2*x*x - 1)) <= 2e-18)
            #expect(abs(try SpecialFunctions.chebyshevT(3, x) - (4*x*x*x - 3*x)) <= 6e-18)

            #expect(try SpecialFunctions.chebyshevU(0, x) == 1)
            #expect(abs(try SpecialFunctions.chebyshevU(1, x) - (2*x)) <= 1e-18)
            #expect(abs(try SpecialFunctions.chebyshevU(2, x) - (4*x*x - 1)) <= 2e-18)
            #expect(abs(try SpecialFunctions.chebyshevU(3, x) - (8*x*x*x - 4*x)) <= 8e-18)
        }
    }

    @Test("Trig identities on |x| ≤ 1 (Float80)")
    func trigIdentitiesFloat80() async throws {
        let thetas: [Float80] = stride(from: Float80(0.2), through: Float80(2.8), by: Float80(0.3)).map { $0 }
        let ns = [0, 1, 2, 3, 5, 7]
        for θ in thetas {
            let x = cos(θ)
            for n in ns {
                let Tn = try SpecialFunctions.chebyshevT(n, x)
                let Un = try SpecialFunctions.chebyshevU(n, x)
                let Tref = cos(Float80(n) * θ)
                let Uref = sin(Float80(n + 1) * θ) / sin(θ)
                #expect(abs(Tn - Tref) <= 1e-18)
                #expect(abs(Un - Uref) <= 2e-18)
            }
        }
    }

    @Test("Endpoint identities (Float80)")
    func endpointsFloat80() async throws {
        for n in 0...16 {
            let T1 = try SpecialFunctions.chebyshevT(n, 1.0 as Float80)
            let Tm1 = try SpecialFunctions.chebyshevT(n, -1.0 as Float80)
            let U1 = try SpecialFunctions.chebyshevU(n, 1.0 as Float80)
            let Um1 = try SpecialFunctions.chebyshevU(n, -1.0 as Float80)
            let s: Float80 = (n % 2 == 0) ? 1 : -1
            #expect(abs(T1 - 1) <= 1e-18)
            #expect(abs(Tm1 - s) <= 1e-18)
            #expect(abs(U1 - Float80(n + 1)) <= 1e-18)
            #expect(abs(Um1 - s * Float80(n + 1)) <= 1e-18)
        }
    }
    #endif

    // MARK: - Recurrence helper

    @Test("chebyshev_next builds T and U sequences (Double)")
    func chebyshevNextDouble() async throws {
        let x: Double = 0.5
        // T: T0=1, T1=x
        let T2 = try SpecialFunctions.chebyshev_next(x, 1.0, x)
        let T3 = try SpecialFunctions.chebyshev_next(T2, x, x)
        #expect(abs(T2 - (-0.5)) <= 1e-15)
        #expect(abs(T3 - (-1.0)) <= 1e-15)

        // U: U0=1, U1=2x
        let U1 = 2.0 * x
        let U2 = try SpecialFunctions.chebyshev_next(U1, 1.0, x)
        let U3 = try SpecialFunctions.chebyshev_next(U2, U1, x)
        #expect(abs(U2 - 0.0) <= 1e-15)
        #expect(abs(U3 - (-1.0)) <= 1e-15)
    }

    @Test("chebyshev_next error paths")
    func chebyshevNextErrors() async throws {
        #expect(throws: SpecialFunctionError.parameterNotFinite(name: "Pn")) {
            _ = try SpecialFunctions.chebyshev_next(Double.nan, 1.0, 0.0 as Double)
        }
        #expect(throws: SpecialFunctionError.parameterNotFinite(name: "Pn1")) {
            _ = try SpecialFunctions.chebyshev_next(1.0 as Double, Double.infinity, 0.0 as Double)
        }
        #expect(throws: SpecialFunctionError.parameterNotFinite(name: "x")) {
            _ = try SpecialFunctions.chebyshev_next(1.0 as Double, 1.0 as Double, -Double.infinity)
        }
        #expect(throws: SpecialFunctionError.parameterNotFinite(name: "Pn")) {
            _ = try SpecialFunctions.chebyshev_next(Float.nan, 1.0 as Float, 0.0 as Float)
        }
        #if arch(x86_64)
        #expect(throws: SpecialFunctionError.parameterNotFinite(name: "x")) {
            _ = try SpecialFunctions.chebyshev_next(1.0 as Float80, 1.0 as Float80, Float80.nan)
        }
        #endif
    }

    // MARK: - Clenshaw recurrence

    @Test("Clenshaw: basic series and half-weight handling (Double)")
    func clenshawDouble() async throws {
        let x: Double = 0.3

        // Default half-weight: S(x) = c0/2 + c1 T1(x) + c2 T2(x) + ...
        // With c = [2, 0, 0], expect constant 1.
        #expect(abs(try SpecialFunctions.chebyshevClenshawRecurrence([2.0], x: x) - 1.0) <= 2e-15)

        // With c = [0, 1], S(x) = T1(x) = x
        #expect(abs(try SpecialFunctions.chebyshevClenshawRecurrence([0.0, 1.0], x: x) - x) <= 2e-15)

        // With c = [0, 0, 1], S(x) = T2(x) = 2x^2 - 1
        let t2 = 2*x*x - 1
        #expect(abs(try SpecialFunctions.chebyshevClenshawRecurrence([0.0, 0.0, 1.0], x: x) - t2) <= 3e-15)

        // No half-weight: S(x) = Σ c_k T_k(x)
        // c = [1, 0, 0] → constant 1
        #expect(abs(try SpecialFunctions.chebyshevClenshawRecurrence([1.0], x: x, halfWeightC0: false) - 1.0) <= 2e-15)
        // c = [0, 1] → T1(x) = x
        #expect(abs(try SpecialFunctions.chebyshevClenshawRecurrence([0.0, 1.0], x: x, halfWeightC0: false) - x) <= 2e-15)
    }

    @Test("Clenshaw: basic series and half-weight handling (Float)")
    func clenshawFloat() async throws {
        let x: Float = 0.3
        #expect(abs(try SpecialFunctions.chebyshevClenshawRecurrence([2.0 as Float], x: x) - 1.0) <= 2e-6)
        #expect(abs(try SpecialFunctions.chebyshevClenshawRecurrence([0.0, 1.0], x: x) - x) <= 2e-6)
        let t2 = 2*x*x - 1
        #expect(abs(try SpecialFunctions.chebyshevClenshawRecurrence([0.0, 0.0, 1.0], x: x) - t2) <= 4e-6)

        #expect(abs(try SpecialFunctions.chebyshevClenshawRecurrence([1.0], x: x, halfWeightC0: false) - 1.0) <= 2e-6)
        #expect(abs(try SpecialFunctions.chebyshevClenshawRecurrence([0.0, 1.0], x: x, halfWeightC0: false) - x) <= 2e-6)
    }

    #if arch(x86_64)
    @Test("Clenshaw: basic series and half-weight handling (Float80)")
    func clenshawFloat80() async throws {
        let x: Float80 = 0.3
        #expect(abs(try SpecialFunctions.chebyshevClenshawRecurrence([2.0 as Float80], x: x) - 1.0) <= 1e-18)
        #expect(abs(try SpecialFunctions.chebyshevClenshawRecurrence([0.0, 1.0], x: x) - x) <= 1e-18)
        let t2 = 2*x*x - 1
        #expect(abs(try SpecialFunctions.chebyshevClenshawRecurrence([0.0, 0.0, 1.0], x: x) - t2) <= 2e-18)

        #expect(abs(try SpecialFunctions.chebyshevClenshawRecurrence([1.0], x: x, halfWeightC0: false) - 1.0) <= 1e-18)
        #expect(abs(try SpecialFunctions.chebyshevClenshawRecurrence([0.0, 1.0], x: x, halfWeightC0: false) - x) <= 1e-18)
    }
    #endif

    @Test("Clenshaw: empty coefficients returns 0 and error paths")
    func clenshawErrors() async throws {
        // Empty coefficients => 0
        #expect(try SpecialFunctions.chebyshevClenshawRecurrence([Double](), x: 0.2) == 0)
        #expect(try SpecialFunctions.chebyshevClenshawRecurrence([Float](), x: 0.2) == 0)
        #if arch(x86_64)
        #expect(try SpecialFunctions.chebyshevClenshawRecurrence([Float80](), x: 0.2) == 0)
        #endif

        // Non-finite x
        #expect(throws: SpecialFunctionError.parameterNotFinite(name: "x")) {
            _ = try SpecialFunctions.chebyshevClenshawRecurrence([1.0], x: Double.nan)
        }
        #expect(throws: SpecialFunctionError.parameterNotFinite(name: "x")) {
            _ = try SpecialFunctions.chebyshevClenshawRecurrence([1.0 as Float], x: Float.infinity)
        }
        #if arch(x86_64)
        #expect(throws: SpecialFunctionError.parameterNotFinite(name: "x")) {
            _ = try SpecialFunctions.chebyshevClenshawRecurrence([1.0 as Float80], x: -Float80.infinity)
        }
        #endif

        // Non-finite coefficients
        #expect(throws: SpecialFunctionError.parameterNotFinite(name: "coefficients")) {
            _ = try SpecialFunctions.chebyshevClenshawRecurrence([1.0, Double.nan], x: 0.2)
        }
        #expect(throws: SpecialFunctionError.parameterNotFinite(name: "coefficients")) {
            _ = try SpecialFunctions.chebyshevClenshawRecurrence([1.0 as Float, -Float.infinity], x: 0.2 as Float)
        }
        #if arch(x86_64)
        #expect(throws: SpecialFunctionError.parameterNotFinite(name: "coefficients")) {
            _ = try SpecialFunctions.chebyshevClenshawRecurrence([1.0 as Float80, Float80.infinity], x: 0.2 as Float80)
        }
        #endif
    }

    // MARK: - Error paths for T and U

    @Test("chebyshevT/chebyshevU error paths")
    func functionErrors() async throws {
        // n < 0
        #expect(throws: SpecialFunctionError.parameterNotPositive(name: "n")) {
            _ = try SpecialFunctions.chebyshevT(-1, 0.3 as Double)
        }
        #expect(throws: SpecialFunctionError.parameterNotPositive(name: "n")) {
            _ = try SpecialFunctions.chebyshevU(-5, 0.3 as Float)
        }
        #if arch(x86_64)
        #expect(throws: SpecialFunctionError.parameterNotPositive(name: "n")) {
            _ = try SpecialFunctions.chebyshevT(-2, 0.3 as Float80)
        }
        #endif

        // Non-finite x
        #expect(throws: SpecialFunctionError.parameterNotFinite(name: "x")) {
            _ = try SpecialFunctions.chebyshevT(2, Double.nan)
        }
        #expect(throws: SpecialFunctionError.parameterNotFinite(name: "x")) {
            _ = try SpecialFunctions.chebyshevU(3, Float.infinity)
        }
        #if arch(x86_64)
        #expect(throws: SpecialFunctionError.parameterNotFinite(name: "x")) {
            _ = try SpecialFunctions.chebyshevU(1, -Float80.infinity)
        }
        #endif
    }
}
