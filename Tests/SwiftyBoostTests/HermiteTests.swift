import Testing
@testable import SwiftyBoost
import CBoostBridge

@Suite("Hermite Polynomials")
struct HermiteTests {

    @Test("Double evaluations match bridge and recurrence")
    func doubleMatchesBridge() throws {
        let xs: [Double] = [-2.0, -0.5, 0.0, 0.75, 1.2]
        for n in 0...8 {
            let nu = UInt32(n)
            for x in xs {
                let expected = bs_hermite_d(nu, x)
                let got = try SpecialFunctions.hermite(n, x)
                #expect(got == expected)
            }
        }
        // Recurrence: H_{n+1} from H_n and H_{n-1}
        for x in xs {
            var hnm1 = try SpecialFunctions.hermite(0, x)
            var hn = try SpecialFunctions.hermite(1, x)
            for n in 1...6 {
                let expectedNext = bs_hermite_next_d(UInt32(n), x, hn, hnm1)
                let gotNext = try SpecialFunctions.hermiteNext(n, x, hn: hn, hnm1: hnm1)
                #expect(gotNext == expectedNext)
                let direct = try SpecialFunctions.hermite(n + 1, x)
                #expect(abs(gotNext - direct) <= 1e-12)
                hnm1 = hn
                hn = gotNext
            }
        }
    }

    @Test("Float overloads mirror bridge")
    func floatMatchesBridge() throws {
        let xs: [Float] = [-1.5, -0.25, 0.0, 0.9, 1.5]
        for n in 0...6 {
            let nu = UInt32(n)
            for x in xs {
                let expected = bs_hermite_f(nu, x)
                let got = try SpecialFunctions.hermite(n, x)
                #expect(got == expected)
            }
        }
        for x in xs {
            var hnm1 = try SpecialFunctions.hermite(0, x)
            var hn = try SpecialFunctions.hermite(1, x)
            for n in 1...5 {
                let expectedNext = bs_hermite_next_f(UInt32(n), x, hn, hnm1)
                let gotNext = try SpecialFunctions.hermiteNext(n, x, hn: hn, hnm1: hnm1)
                #expect(gotNext == expectedNext)
                let direct = try SpecialFunctions.hermite(n + 1, x)
                #expect(abs(Double(gotNext) - Double(direct)) <= 1e-5)
                hnm1 = hn
                hn = gotNext
            }
        }
    }

    #if arch(x86_64) || arch(i386)
    @Test("Float80 overloads mirror bridge")
    func float80MatchesBridge() throws {
        let xs: [Float80] = [-1.25, -0.5, 0.0, 0.6, 1.1]
        for n in 0...5 {
            let nu = UInt32(n)
            for x in xs {
                let expected = bs_hermite_l(nu, x)
                let got = try SpecialFunctions.hermite(n, x)
                #expect(got == expected)
            }
        }
        for x in xs {
            var hnm1 = try SpecialFunctions.hermite(0, x)
            var hn = try SpecialFunctions.hermite(1, x)
            for n in 1...4 {
                let expectedNext = bs_hermite_next_l(UInt32(n), x, hn, hnm1)
                let gotNext = try SpecialFunctions.hermiteNext(n, x, hn: hn, hnm1: hnm1)
                #expect(gotNext == expectedNext)
                let direct = try SpecialFunctions.hermite(n + 1, x)
                #expect(abs(gotNext - direct) <= 1e-15)
                hnm1 = hn
                hn = gotNext
            }
        }
    }
    #endif

    @Test("Rejects invalid parameters")
    func rejectsInvalidInputs() {
        #expect(throws: SpecialFunctionError<Double>.parameterNotPositive(name: "n", value: -1)) {
            _ = try SpecialFunctions.hermite(-1, 0.2 as Double)
        }
        do {
            _ = try SpecialFunctions.hermite(2, Double.nan)
            #expect(Bool(false), "Expected parameterNotFinite for x=nan")
        } catch let error as SpecialFunctionError<Double> {
            if case let .parameterNotFinite(name: actualName, value: value) = error {
                #expect(actualName == "x")
                #expect(value.isNaN)
            } else {
                #expect(Bool(false), "Unexpected error: \(error)")
            }
        } catch {
            #expect(Bool(false), "Unexpected error type: \(error)")
        }
        #expect(throws: SpecialFunctionError<Double>.parameterNotPositive(name: "n", value: 0)) {
            _ = try SpecialFunctions.hermiteNext(0, 0.2 as Double, hn: 1.0, hnm1: 1.0)
        }
        do {
            _ = try SpecialFunctions.hermiteNext(1, Double.infinity, hn: 2.0, hnm1: 1.0)
            #expect(Bool(false), "Expected parameterNotFinite for x=∞")
        } catch let error as SpecialFunctionError<Double> {
            if case let .parameterNotFinite(name: actual, value: value) = error {
                #expect(actual == "x")
                #expect(value.isInfinite)
            } else {
                #expect(Bool(false), "Unexpected error: \(error)")
            }
        } catch {
            #expect(Bool(false), "Unexpected error type: \(error)")
        }
        do {
            _ = try SpecialFunctions.hermiteNext(2, 0.5 as Double, hn: Double.nan, hnm1: 1.0)
            #expect(Bool(false), "Expected parameterNotFinite for hn=nan")
        } catch let error as SpecialFunctionError<Double> {
            if case let .parameterNotFinite(name: actual, value: value) = error {
                #expect(actual == "hn")
                #expect(value.isNaN)
            } else {
                #expect(Bool(false), "Unexpected error: \(error)")
            }
        } catch {
            #expect(Bool(false), "Unexpected error type: \(error)")
        }
        do {
            _ = try SpecialFunctions.hermiteNext(2, 0.5 as Double, hn: 1.0, hnm1: Double.infinity)
            #expect(Bool(false), "Expected parameterNotFinite for hnm1=∞")
        } catch let error as SpecialFunctionError<Double> {
            if case let .parameterNotFinite(name: actual, value: value) = error {
                #expect(actual == "hnm1")
                #expect(value.isInfinite)
            } else {
                #expect(Bool(false), "Unexpected error: \(error)")
            }
        } catch {
            #expect(Bool(false), "Unexpected error type: \(error)")
        }
    }
}
