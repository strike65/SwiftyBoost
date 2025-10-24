import Testing
@testable import SwiftyBoost
import CBoostBridge

@Suite("Jacobi Zeta Function")
struct JacobiZetaTests {

    @Test("Double results match bridge")
    func doubleMatchesBridge() throws {
        let phis: [Double] = [-1.2, -0.4, 0.0, 0.5, 1.1]
        let ks: [Double] = [0.0, 0.2, 0.75, 0.99]
        for phi in phis {
            for k in ks {
                let expected = bs_jacobi_zeta_d(k, phi)
                let got = try SpecialFunctions.jacobiZeta(phi, modulus: k)
                #expect(got == expected)
            }
        }
    }

    @Test("Float results match bridge")
    func floatMatchesBridge() throws {
        let phis: [Float] = [-0.9, -0.25, 0.0, 0.45, 1.0]
        let ks: [Float] = [0.1, 0.5, 0.95]
        for phi in phis {
            for k in ks {
                let expected = bs_jacobi_zeta_f(k, phi)
                let got = try SpecialFunctions.jacobiZeta(phi, modulus: k)
                #expect(got == expected)
            }
        }
    }

    #if arch(x86_64) || arch(i386)
    @Test("Float80 results match bridge")
    func float80MatchesBridge() throws {
        let phis: [Float80] = [-0.7, -0.3, 0.0, 0.6, 1.2]
        let ks: [Float80] = [0.0, 0.6, 0.99]
        for phi in phis {
            for k in ks {
                let expected = bs_jacobi_zeta_l(k, phi)
                let got = try SpecialFunctions.jacobiZeta(phi, modulus: k)
                #expect(got == expected)
            }
        }
    }
    #endif

    @Test("Rejects non-finite inputs")
    func rejectsInvalid() {
        #expect(throws: SpecialFunctionError<Double>.parameterNotFinite(name: "phi", value: .infinity)) {
            _ = try SpecialFunctions.jacobiZeta(.infinity, modulus: 0.5)
        }
        do {
            _ = try SpecialFunctions.jacobiZeta(0.5, modulus: .nan)
            #expect(Bool(false), "Expected parameterNotFinite for k")
        } catch let error as SpecialFunctionError<Double> {
            if case let .parameterNotFinite(name: actual, value: value) = error {
                #expect(actual == "k")
                #expect(value.isNaN)
            } else {
                #expect(Bool(false), "Unexpected error: \(error)")
            }
        } catch {
            #expect(Bool(false), "Unexpected error type: \(error)")
        }
    }
}
